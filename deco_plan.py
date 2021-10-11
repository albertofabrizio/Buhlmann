import argparse
import numpy as np
from tabulate import tabulate
from tools import Buhlmann as balgo

########################## Parsing user defined input ##########################
parser = argparse.ArgumentParser(description='This program computes the decompression obligations based on the Buhlmann algorithm and IANTD standards.')

#parser.add_argument('--rep', action='store_true', help="Wheter this is a repetitive dive.")
parser.add_argument('--altitude', type=float, dest='h_on_sea', default=0.0, help="Altitude on sea-level for the dive [m] [default: 0.0]")
parser.add_argument('--alveolar', type=str, dest='alveolar', default='buhl', help="Which water vapor pressure in the lungs to use. ['buhl', 'schrein', 'navy']")
parser.add_argument('--depth', type=float, dest='tdepth', required=True, help="Target maximum depth [m]")
parser.add_argument('--time', type=float, dest='runT', required=True, help="Run time [min] at which you desire to quit the target maximum depth.")
parser.add_argument('--fo2', type=float, dest='fo2', required=True, nargs='+', help="Fractions of oxygen in breathing gas (one value per tank).")
parser.add_argument('--fhe', type=float, dest='fhe', required=True, nargs='+', help="Fractions of helium in breathing gas (one value per tank).")
parser.add_argument('--SAC', type=float, dest='sac', default='20', help="Surface air consumption [l/min] [default: 20].")
parser.add_argument('--glow', type=float, dest='gf_low', default='0.75', help="Gradient factor (Low) [default: 0.75].")
parser.add_argument('--ghigh', type=float, dest='gf_hi', default='0.75', help="Gradient factor (High) [default: 0.75].")
parser.add_argument('--last', type=float, dest='last_deco', default='6', help="Last deco stop [m] [default: 6].")
parser.add_argument('--debug', action='store_true', help="Print Debug Info.")


args = parser.parse_args()

########################## Helper Functions ######################## 

def convert_to_depth_Abs(pressure):
    # Initialize environment conditions.
    env = Environment(args.h_on_sea)

    # Compute ambient pressure given altitude on the sea level.
    pAlt = env.altitude()

    depth = (float(pressure) - pAlt) * 10.0 
    return depth

def convert_to_bar_Abs(depth):
    # Initialize environment conditions.
    env = Environment(args.h_on_sea)

    # Compute ambient pressure given altitude on the sea level.
    pAlt = env.altitude()

    bar = pAlt + float(depth) / 10.0
    return bar

########################## Main program and Objects ######################## 


class Constants():
    """ Phisiological, chemical, physical constants."""

    # For more info about Alveolar WVP see https://www.shearwater.com/wp-content/uploads/2012/08/Introductory-Deco-Lessons.pdf
    AlveolarWVP_Buhlmann = 0.0627 # Alveolar water vapor partial pressure at surface. Original from Buhlmann paper (assuming respiratory quotient of 1.0).
    AlveolarWVP_Schreiner = 0.0493 # Alveolar water vapor partial pressure at surface. From Schreiner (assuming respiratory quotient of 0.8).
    AlveolarWVP_Navy = 0.0567 # Alveolar water vapor partial pressure at surface. From U.S. Navy (assuming respiratory quotient of 0.9)

    if (args.alveolar == 'buhl'):
        AlveolarWVP = AlveolarWVP_Buhlmann

    elif (args.alveolar == 'schrein'):
        AlveolarWVP = AlveolarWVP_Schreiner

    else:
        AlveolarWVP = AlveolarWVP_Navy
    
    
    surfacePressure = 1.01325  # Pressure at sea level [bar].
    surfaceTemperature = 20 # Temperature at surface [C].

    AirHelium = 0.0  # Fraction of He in standard air.
    AirNitrogen = 0.7808  # Fraction of N2 in standard air.
    AirOxygen = 0.2095  # Fraction of O2 in standard air.
    AirArgon = 0.00934  # Fraction of Ar in standard air.

    AirInertGasFraction = AirHelium + AirNitrogen + AirArgon


    _M = 0.02896 # Molar mass of earth's air
    _g = 9.807 # Earth's gravitational acceleration
    _R = 8.3143 # Ideal gas constant
    _T = 288.15 # Standard temperature


class Environment():

    def __init__(self, h_on_sea):
        self.h_on_sea = h_on_sea

    def altitude(self):
        """ Returns the total pressure of air given the height over sea-level in meters."""
        
        # Constants for the barometric formula
        cnst = (Constants._M*Constants._g)/(Constants._R*Constants._T)

        return Constants.surfacePressure * np.exp( -cnst * self.h_on_sea)

class Compartments():
    ''' Object containing the Buhlmann compartments. '''

    def __init__(self, params):

        self.d = Constants.surfacePressure 
        self.speed_descent = 20 / 10.0
        self.speed_deep = 9 / 10.0
        self.speed_shallow = 3 / 10.0

        self.compartments = np.zeros((16,3)) # Column are pN2, pHe and pInert

        assert len(args.fo2) == len(args.fhe), "The list of Oxygen and Helium fractions does not have the same length."
        self.fo2 = np.array(args.fo2)
        self.fhe = np.array(args.fhe)
        self.fn2 = np.ones_like(self.fo2) - self.fo2 - self.fhe

        self.mod = self.mod_o2()    
        self.current_gas = 0
 
        self.ht_n2 = params[:,0]
        self.ht_he = params[:,1]

        self.a_n2 = params[:,2] #* Constants.surfacePressure
        self.b_n2 = params[:,3]
        self.a_he = params[:,4] #* Constants.surfacePressure
        self.b_he = params[:,5]

        self.ascent_ceil = np.zeros(16)
        self.deco_stop = 0
        self.deco_stops = []
        self.max_loading = 0
        self.sac = args.sac
        self.volume = []

        self.GF = 1.0
        self.first_stop = 0.0
        self.first_stop_flag = True

        self.deco_profile = {}


    def initialize(self):
        """ Initialize inert gasses partial pressures and total inert gas pressure for all the 16 tissues."""

        # Communicate with user.
        print("Initializing Tissue Compartments. Altitude in [m]:", args.h_on_sea)
        
        # Initialize environment conditions.
        env = Environment(args.h_on_sea)

        # Compute ambient pressure given altitude on the sea level.
        self.palt = env.altitude()

        # Evaluate N2 partial pressure and initialize the compartments.
        pN2 = (self.palt - Constants.AlveolarWVP) * Constants.AirNitrogen
        pTot = pN2 + Constants.AirHelium
    
        for comp_idx in range(16):
            self.compartments[comp_idx,0] = pN2
            self.compartments[comp_idx,1] = Constants.AirHelium
            self.compartments[comp_idx,2] = pTot

        if args.debug:
            self.print_comp()

        # Check if there is a gas appropriate for target depth and use it first. (Assumption: no hypoxic)
        self.current_gas = np.argmax(self.mod)
        assert (self.mod[self.current_gas] - self.palt) * 10.0 >= args.tdepth, "Planned depth exceeds bottom gas MOD."
        
        # Go to target depth
        self.constant_speed(0, args.tdepth)

    def mod_o2(self):
        return np.divide(1.6,self.fo2)

    def isobaric_CD_ratio(self, oldN2, newN2, oldHe, newHe):

        test = (newHe - oldHe) * 100 / 5

        if np.round(np.abs(test),1) < np.round((newN2 - oldN2) * 100, 1) and np.round((newN2 - oldN2),1) > 0 and np.round(test,1) < 0:
            print("Warning: the increase in N2 content exceeds 5 times the decrease in He. Possible ICD warning.")
            print("N2 Increase:", np.round((newN2 - oldN2) * 100,1),"%")
            print("He Decrease:", np.round((newHe - oldHe) * 100,1),"%")

    def check_better_gas(self):

        i = 0
        tmp_old = 2 # If tmp is bigger than 1.61325, the argument of the if statement will be anyway False
        for fo2 in self.fo2:
            tmp = 1.61325 - self.d * fo2
            if tmp >= 0 and tmp_old > tmp:
                tmp_old = tmp
                better = i                
            i+=1

        if better != self.current_gas:
            print("Switching gas at a depth of:", np.round((self.d - self.palt)*10.0), "m for fO2 of: ", self.fo2[better])

            # If there is Helium in the old mix check for ICD ratio upon switch.
            if self.fhe[self.current_gas] != 0:
                self.isobaric_CD_ratio(self.fn2[self.current_gas], self.fn2[better], self.fhe[self.current_gas], self.fhe[better])

            self.current_gas = better

    def constant_depth(self, depth, time):

        self.d = convert_to_bar_Abs(depth)

        # Check if a better gas exist
        self.check_better_gas()

        for comp_idx in range(16):

            p0_n2 = self.compartments[comp_idx,0]
            p0_he = self.compartments[comp_idx,1]

            pi_n2 = (self.d - Constants.AlveolarWVP) * self.fn2[self.current_gas]
            pN2 = p0_n2 + (pi_n2 - p0_n2)*(1-2.0**(-time/self.ht_n2[comp_idx]))

            pi_he = (self.d - Constants.AlveolarWVP) * self.fhe[self.current_gas]
            pHe = p0_he + (pi_he - p0_he)*(1-2.0**(-time/self.ht_he[comp_idx]))

            self.compartments[comp_idx,0] = pN2
            self.compartments[comp_idx,1] = pHe
            self.compartments[comp_idx,2] = pN2 + pHe

        if args.debug:
            self.print_comp()

        # Once we know the saturation of the tissues, check the ascent ceiling
        self.check_ascent_ceiling()

    def constant_speed(self, depth_i, depth_f):

        d1 = convert_to_bar_Abs(depth_i)
        d2 = convert_to_bar_Abs(depth_f)

        for comp_idx in range(16):
            p0_n2 = self.compartments[comp_idx,0]
            p0_he = self.compartments[comp_idx,1]

            pi0_n2 = (d1 - Constants.AlveolarWVP) * self.fn2[self.current_gas]
            pi0_he = (d1 - Constants.AlveolarWVP) * self.fhe[self.current_gas]

            if (depth_f > depth_i):
                RN2a = self.speed_descent * self.fn2[self.current_gas]
                RHe = self.speed_descent * self.fhe[self.current_gas]
                t = (d2-d1) / self.speed_descent

            else: 
                if (depth_f - depth_i < -3):
                    RN2a = -self.speed_deep * self.fn2[self.current_gas]
                    RHe = -self.speed_deep * self.fhe[self.current_gas]
                    t = (d2-d1) / -self.speed_deep
                else:
                    RN2a = -self.speed_shallow * self.fn2[self.current_gas]
                    RHe = -self.speed_shallow * self.fhe[self.current_gas]
                    t = (d2-d1) / -self.speed_shallow
                    
            kN2 = np.log(2) / self.ht_n2[comp_idx]
            kHe = np.log(2) / self.ht_he[comp_idx]

            pN2 = pi0_n2 + RN2a*(t - (1/kN2)) - (pi0_n2 - p0_n2 - (RN2a / kN2))* np.exp(-kN2*t)
            pHe = pi0_he + RHe*(t - (1/kHe)) - (pi0_he - p0_he - (RHe / kHe))* np.exp(-kHe*t)

            self.compartments[comp_idx,0] = pN2
            self.compartments[comp_idx,1] = pHe
            self.compartments[comp_idx,2] = pN2 + pHe

        if args.debug:
            self.print_comp()


    def compute_GF(self):

        ''' 
        The idea is that we draw a line that connects GF_high at the surface to GF_low at the first stop. 
        The equation of this line is GF(d) = (d - pAlt)/(pFirst - pAlt) * (GF_low - GF_high) + GF_high.
        Of course this number only make sense when d <= pFirst...
        '''

        self.GF = (self.d - self.palt)/(self.first_stop - self.palt) * (args.gf_low - args.gf_hi) + args.gf_hi
        

    def compute_stop_in_m(self):
        """ Converts the inert gas loading to ceiling stop and terminates the program if NDL is not reached. """

        # Check among all tissues, which one has the highest loading.
        self.max_loading = np.max(self.ascent_ceil)
                
        # Convert the highest loading to a depth increment of 3 m.
        self.deco_stop = convert_to_depth_Abs(self.max_loading)
        self.deco_stop = 3 * np.ceil(self.deco_stop/3)
                
        if self.deco_stop < 2:
            self.deco_stop = 0.0

        self.deco_stops.append(self.deco_stop)

        # Kill the program if there is no ceiling for the dive.
        if (len(self.deco_stops) == 1 and self.deco_stops[0] == 0.0):
            print("This is a dive within the NDL. Stop at:")
            print("9 m --- 1 min")
            print("6 m --- 2 min")
            print("3 m --- 3 min")
            exit()


    def check_ascent_ceiling(self):
        """ Computes the ceiling depth given the current inert gas loading. """
        
        # If it is not the first time that we compute the ceiling... Then we need to know the GF for the current depth.
        if not self.first_stop_flag:
            self.compute_GF()

        # For each compartment
        for comp_idx in range(16):

            # Compute Buhlmann M-value line parameters A and B with both N2 and He
            A = ((self.a_n2[comp_idx] * self.compartments[comp_idx,0]) + (self.a_he[comp_idx] * self.compartments[comp_idx,1]))/(self.compartments[comp_idx,2])
            B = ((self.b_n2[comp_idx] * self.compartments[comp_idx,0]) + (self.b_he[comp_idx] * self.compartments[comp_idx,1]))/(self.compartments[comp_idx,2])

            # If it is the first time that we compute the ceiling depth: use GF_low.
            if self.first_stop_flag:
                self.ascent_ceil[comp_idx] = ((self.compartments[comp_idx,0]+self.compartments[comp_idx,1]) - args.gf_low * A) / ( args.gf_low/B - args.gf_low + 1)
            
            # Otherwise use self.GF computed precendently.
            else:
                self.ascent_ceil[comp_idx] = ((self.compartments[comp_idx,0]+self.compartments[comp_idx,1]) - self.GF * A) / ( self.GF/B - self.GF + 1)

        # Convert pressures in meters and stop the program if the dive in within the NDL.
        self.compute_stop_in_m()

        # We want our last deco to be at the user-selected depth.
        if (self.deco_stop < args.last_deco and self.deco_stop > 0):
            self.deco_stop = args.last_deco        
        
        # If this was the first time that we computed ceil... we need to attribute the first_stop variable to draw the GF(depth) line.
        if (self.first_stop_flag):    
            self.first_stop_flag = False
            self.first_stop = convert_to_bar_Abs(self.deco_stop)

            # Also we need to initialize the dictionary that will contain all our deco stops.
            for depth in range(0,int(args.tdepth)+1,3):
                self.deco_profile[depth] = 0   

        # First we need to know at which depth we are in [m]
        current = round(convert_to_depth_Abs(self.d))

        # If the current depth is not a deco stop... we need first to go there.
        if current != int(self.deco_stop):

            # First check if you need to stop to switch gas
            stop_press = convert_to_bar_Abs(self.deco_stop)

            i = 0 
            deco_gas_idx = self.current_gas
            for gas in self.mod:
                diff = gas - stop_press
                if diff > 0:
                    deco_gas_idx = i                    
                i += 1

            # If it is the case
            if deco_gas_idx != self.current_gas:
                # Go to MOD of the gas and stay there 2 minutes

                self.deco_stop = round(convert_to_depth_Abs(self.mod[deco_gas_idx])/3)*3

                self.deco_profile[int(self.deco_stop)] += 2
                self.constant_speed(current,self.deco_stop)
                self.check_better_gas()
                # print(self.current_gas)
                self.constant_depth(self.deco_stop, 2)
                return

            else:
                # We are already using the best gas and just go to the next ceiling
                self.constant_speed(current,self.deco_stop)

        # In all cases (whether you moved or not)... you need to stay at the ceiling depth
        if self.deco_stop > 0 :
            # We stop there for a minute and update our compartments
            self.deco_profile[int(self.deco_stop)] += 1
            self.constant_depth(self.deco_stop, 1)        
            return

    def compute_needed_volume(self):
        
        # Add the volume needed for deco
        for key, value in self.deco_profile.items():
            if key > 0:
                if value > 0 :
                    self.volume.append(self.sac * convert_to_bar_Abs(key) * value)

        # Compute the gas needed for descent
        time_descent = np.abs(convert_to_bar_Abs(args.tdepth)- convert_to_bar_Abs(0)) / self.speed_descent
        desc_gas = self.sac * (convert_to_bar_Abs(args.tdepth) + convert_to_bar_Abs(0)) / 2 * time_descent

        # Compute the gas needed for bottom time at max depth
        bottom_gas = self.sac * convert_to_bar_Abs(args.tdepth) * args.runT
        self.volume.append(desc_gas + bottom_gas)
            

    def print_comp(self):
        """ Print the compartment saturation in a beautiful table. """

        print("")
        print(tabulate(self.compartments, headers=['pN2 [bar]', 'pHe [bar]', 'p_Inert [bar]']))

    def print_profile(self):
        """ Print the computed dive profile in a beautiful table. """
        
        tmp = []
        
        # Enforce IANTD standards.
        if self.deco_profile[9] < 1:
            self.deco_profile[9] = 1

        if self.deco_profile[6] < 2:
            self.deco_profile[6] = 2

        if self.deco_profile[3] < 3:
            self.deco_profile[3] = 3

        self.compute_needed_volume()

        i=0
        for key, value in self.deco_profile.items():
            if key > 0:
                if value > 0 :
                    tmp.append([key, value, "{0:.2f}".format(self.volume[i])])
                    i+=1

        tmp.append([args.tdepth, args.runT, "{0:.2f}".format(self.volume[-1],3)])
        

        print("")
        print(tabulate(tmp, headers=['Depth [m]', 'Stop [min]', 'Volume Gas [l]']))

        


########################## Main ########################

def main():

    # Load the parameters
    params = balgo.ZHL16c_1b()
    
    # Initialize all the tissues to air partial pressures.
    tissues = Compartments(params)
    tissues.initialize()

    # Stay at target depth for given bottom time
    tissues.constant_depth(args.tdepth, args.runT)

    # Print the deco profile
    tissues.print_profile()


if __name__ == "__main__":
    main()
