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
parser.add_argument('--speed_d', type=float, dest='speed_descent', default='20', help="Speed of descent [m/min] [default: 18].")
parser.add_argument('--speed_a_deep', type=float, dest='speed_deep', default='9', help="Speed of ascent from depth to first stop [m/min] [default: 9].")
parser.add_argument('--speed_a_shallow', type=float, dest='speed_shallow', default='3', help="Speed of ascent from stop to stop [m/min] [default: 3].")
parser.add_argument('--mod', type=list, nargs=2, dest='mod_user', default=[1.4,1.6], help="User defined ppO2 at maximum operating depth and at deco [bar] [default: 1.4, 1.6].")


parser.add_argument('--debug', action='store_true', help="Print Debug Info.")


args = parser.parse_args()

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
        """Initialize the altitude on sea level"""
        self.h_on_sea = h_on_sea

    def altitude(self):
        """ Returns the total pressure of air given the height over sea-level in meters."""
        
        # Constants for the barometric formula
        cnst = (Constants._M*Constants._g)/(Constants._R*Constants._T)

        return Constants.surfacePressure * np.exp( -cnst * self.h_on_sea)

class Compartments():
    ''' Object containing the Buhlmann compartments. '''

    def __init__(self, params):

        # Initialize environment
        env = Environment(args.h_on_sea)
        self.palt = env.altitude()

        # Set initial pressure (depth) to surface absolute pressure
        self.d = self.palt 

        # Set ascent and descent speeds in [bar]
        self.speed_descent = args.speed_descent / 10.0
        self.speed_deep = -1 * args.speed_deep / 10.0
        self.speed_shallow = -1 * args.speed_shallow / 10.0

        # Initialize compartments array
        self.compartments = np.zeros((16,3)) # Column are pN2, pHe and pInert

        # Initialize gases arrays and check for sanity
        assert len(args.fo2) == len(args.fhe), "The list of Oxygen and Helium fractions does not have the same length."
        self.fo2 = np.array(args.fo2)
        self.fhe = np.array(args.fhe)
        self.fn2 = np.ones_like(self.fo2) - self.fo2 - self.fhe
        self.gas_labels = [str(int(self.fo2[i]*100)).zfill(2)+"/"+str(int(self.fhe[i]*100)).zfill(2) for i in range(len(self.fo2))]

        # Compute the MOD and set global gas variables
        self.mod_depth = self.mod_o2(isdeco = False)    
        self.mod_deco = self.mod_o2(isdeco = True)    
        self.current_gas = 0
 
        # Load Buhlmann half-times and coeffs
        self.ht_n2 = params[:,0]
        self.ht_he = params[:,1]

        self.a_n2 = params[:,2] 
        self.b_n2 = params[:,3]
        self.a_he = params[:,4] 
        self.b_he = params[:,5]

        # Set SAC and gas volume list
        self.sac = args.sac
        self.gas_consumed_ascent_deep = 0

        # Set global ceiling arrays
        self.ascent_ceil_bar = np.zeros(16)
        self.tol = np.zeros(16)
        self.current_deco_stop = 0
        self.max_loading = 0

        # Time
        self.run = 0
        self.start = 0
        self.start_depth = 0

        # Gradient factors        
        self.GF = 1.0
        self.first_stop = 0.0

        self.dive_profile = {"D":[], "T":[], "R":[], "S":[], "G":[], "V":[]}

########################## Helper Functions ######################## 

    def convert_to_depth(self, pressure):
        """Converts absolute pressure to depth in [m]"""
        
        return (float(pressure) - self.palt) * 10.0 

    def convert_to_press_abs(self, depth):
        """Converts depth in [m] to absolute pressure [bar]"""

        return self.palt + float(depth) / 10.0

    def mod_o2(self, isdeco = False):
        """Returns the maximum operating depth"""

        if isdeco:
            return np.divide(args.mod_user[1],self.fo2)

        else:
            return np.divide(args.mod_user[0],self.fo2)

    def compute_GF(self):
        ''' 
        The idea is that we draw a line that connects GF_high at the surface to GF_low at the first stop. 
        The equation of this line is GF(d) = (d - pAlt)/(pFirst - pAlt) * (GF_low - GF_high) + GF_high.
        Of course this number only make sense when d <= pFirst...
        '''

        last_stop = self.convert_to_press_abs(args.last_deco)
        if self.first_stop != last_stop:
            self.GF = (self.d - last_stop)/(self.first_stop - last_stop) * (args.gf_low - args.gf_hi) + args.gf_hi 
        else:
            self.GF = args.gf_hi


    def print_comp(self):
        """ Print the compartment saturation in a beautiful table. """

        print("")
        print(tabulate(self.compartments, headers=['pN2 [bar]', 'pHe [bar]', 'p_Inert [bar]']))

    def print_profile(self):
        print("")
        print(tabulate(self.dive_profile, headers=['Depth [m]', 'Time [min]', 'Run [min]', 'Stop [min]', 'Gas [O2/He]', 'Volume Gas [l]']))

        TTS = self.dive_profile["R"][-1] - self.dive_profile["R"][1]
        print("")
        print("TTS [min]:", TTS)

        print("")
        print("Total Volume of gas Used:")
        for gas_l in self.gas_labels:
            tot_gas_used = 0 
            for idx, gas in enumerate(self.dive_profile["G"]):
                tot_gas_used += int(gas==gas_l) * self.dive_profile["V"][idx]

            print(gas_l + " [l]:", round(tot_gas_used,2))

    def go_to_switch(self, switch_depth):
        d1 = self.d
        d2 = self.convert_to_press_abs(switch_depth)

        speed = self.speed_deep

        for comp_idx in range(16):
            p0_n2 = self.compartments[comp_idx,0]
            p0_he = self.compartments[comp_idx,1]

            pi0_n2 = (d1 - Constants.AlveolarWVP) * self.fn2[self.current_gas]
            pi0_he = (d1 - Constants.AlveolarWVP) * self.fhe[self.current_gas]

            RN2a = speed * self.fn2[self.current_gas]
            RHe = speed * self.fhe[self.current_gas]
            # The np.round is my choice ... I could leave it as a fraction but I prefer to approximate to the next integer.
            t = np.round((d2-d1) / speed)
                    
            kN2 = np.log(2) / self.ht_n2[comp_idx]
            kHe = np.log(2) / self.ht_he[comp_idx]

            pN2 = pi0_n2 + RN2a*(t - (1/kN2)) - (pi0_n2 - p0_n2 - (RN2a / kN2))* np.exp(-kN2*t)
            pHe = pi0_he + RHe*(t - (1/kHe)) - (pi0_he - p0_he - (RHe / kHe))* np.exp(-kHe*t)

            self.compartments[comp_idx,0] = pN2
            self.compartments[comp_idx,1] = pHe
            self.compartments[comp_idx,2] = pN2 + pHe

        # Add descent to runtime
        self.run += t
        
        # Compute gas needed for descent
        self.gas_consumed_ascent_deep += self.sac * ((d1+d2)/2)/ self.palt * t
        
        # Print Debug info

        if args.debug:
            print("")
            print("****************** Go to Switch ******************")
            print("Starting depth [bar]:", d1)
            print("Final depth [bar]:", d2)
            print("Ascent time [min]:", t)
            print("Ascent speed [bar/min]:", speed)
            print("")
            print("Current Runtime [min]:", self.run)
            print("Current Gas [O2/He]:", self.gas_labels[self.current_gas])
            print("Gas Volume used [l]:", self.gas_consumed_ascent_deep)
            print("Saturation of Tissue Compartments:")
            self.print_comp()
            print("****************** Go to Switch ******************")

        # Check if d2 is still a ceiling (check off-gassing during ascent)      
        self.d = d2

        # Update dive profile
        self.dive_profile["D"].append(str(round(self.convert_to_depth(self.start_depth)))+" -> "+str(round(self.convert_to_depth(d2))))
        self.dive_profile["T"].append(self.start)
        self.dive_profile["R"].append(self.run)
        self.dive_profile["S"].append("-")
        self.dive_profile["G"].append(self.gas_labels[self.current_gas])
        self.dive_profile["V"].append(self.gas_consumed_ascent_deep)

        self.start = self.run
        self.start_depth = round(self.convert_to_depth(d2))

    def isobaric_CD_ratio(self, oldN2, newN2, oldHe, newHe):

        test = (newHe - oldHe) * 100 / 5

        if np.round(np.abs(test),1) < np.round((newN2 - oldN2) * 100, 1) and np.round((newN2 - oldN2),1) > 0 and np.round(test,1) < 0:
            print("Warning: the increase in N2 content exceeds 5 times the decrease in He. Possible ICD warning.")
            print("N2 Increase:", np.round((newN2 - oldN2) * 100,1),"%")
            print("He Decrease:", np.round((newHe - oldHe) * 100,1),"%")

    def stay_at_switch(self):

        time = 2 # Stay 2 minutes at gas switch

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

        # Add descent to runtime
        self.run += time

        # Compute gas needed for descent
        gas_consumed = self.sac * self.d/ self.palt * time

        # Print Debug info
        if args.debug:
            print("")
            print("****************** Bottom ******************")
            print("Current depth [bar]:", self.d)
            print("Time at depth [min]:", time)
            print("")
            print("Current Runtime [min]:", self.run)
            print("Current Gas [O2/He]:", self.gas_labels[self.current_gas])
            print("Gas Volume used [l]:", gas_consumed)
            print("Saturation of Tissue Compartments:")
            self.print_comp()
            print("****************** Bottom ******************")

        # Update Dive Profile
        self.dive_profile["D"].append(str(round(self.convert_to_depth(self.d))))
        self.dive_profile["T"].append(self.start)
        self.dive_profile["R"].append(self.run)
        self.dive_profile["S"].append(int(self.run - self.start))
        self.dive_profile["G"].append(self.gas_labels[self.current_gas])
        self.dive_profile["V"].append(gas_consumed)

    def gas_switch(self, switch_depth):

        # Go to the switch depth with old gas
        self.go_to_switch(switch_depth)

        # Switch Gases
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

        # Stay 2 min
        self.stay_at_switch()

        # go to ceiling loop
        self.start = self.run
        self.start_depth = self.d
        self.check_ascent_ceiling()
        self.ascent_deep(self.current_deco_stop)

########################## Main Deco Algorithm ######################## 

    def initialize(self):
        """ Initialize inert gasses partial pressures and total inert gas pressure for all the 16 tissues."""

        # Communicate with user.
        print("Initializing Tissue Compartments. Altitude in [m]:", args.h_on_sea)
        print("Available Gas Mixes [O2/He]:", self.gas_labels)
        
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

        # Search for appropriate gas mix. (Assumption: no hypoxic)
        self.current_gas = np.argmax(self.mod_depth)
        assert (self.mod_depth[self.current_gas] - self.palt) * 10.0 >= args.tdepth, "Planned depth exceeds bottom gas MOD."

        if args.debug:
            print("")
            print("****************** Initialize ******************")
            print("Surface Pressure [bar]:", self.palt)
            print("Selected Gas for Dive [O2/He]:", self.gas_labels[self.current_gas])
            print("Starting Saturation of Tissue Compartments:")
            self.print_comp()
            print("****************** Initialize ******************")

        self.descent(0, args.tdepth)

    def descent(self, depth_i, depth_f):
        """Compute tissues' inert gas loading during descent -> Schreiner Equation."""

        d1 = self.convert_to_press_abs(depth_i)
        d2 = self.convert_to_press_abs(depth_f)

        for comp_idx in range(16):
            p0_n2 = self.compartments[comp_idx,0]
            p0_he = self.compartments[comp_idx,1]

            pi0_n2 = (d1 - Constants.AlveolarWVP) * self.fn2[self.current_gas]
            pi0_he = (d1 - Constants.AlveolarWVP) * self.fhe[self.current_gas]

            RN2a = self.speed_descent * self.fn2[self.current_gas]
            RHe = self.speed_descent * self.fhe[self.current_gas]

            # The np.round is my choice ... I could leave it as a fraction but I prefer to approximate to the next integer.
            t = np.round((d2-d1) / self.speed_descent)
                    
            kN2 = np.log(2) / self.ht_n2[comp_idx]
            kHe = np.log(2) / self.ht_he[comp_idx]

            pN2 = pi0_n2 + RN2a*(t - (1/kN2)) - (pi0_n2 - p0_n2 - (RN2a / kN2))* np.exp(-kN2*t)
            pHe = pi0_he + RHe*(t - (1/kHe)) - (pi0_he - p0_he - (RHe / kHe))* np.exp(-kHe*t)

            self.compartments[comp_idx,0] = pN2
            self.compartments[comp_idx,1] = pHe
            self.compartments[comp_idx,2] = pN2 + pHe

        # Add descent to runtime
        self.run += t

        # Compute gas needed for descent
        gas_consumed = self.sac * ((d1+d2)/2)/ self.palt * t

        # Print Debug info

        if args.debug:
            print("")
            print("****************** Descent ******************")
            print("Starting depth [bar]:", d1)
            print("Final depth [bar]:", d2)
            print("Descent time [min]:", t)
            print("Descent speed [bar/min]:", self.speed_descent)
            print("")
            print("Current Runtime [min]:", self.run)
            print("Current Gas [O2/He]:", self.gas_labels[self.current_gas])
            print("Gas Volume used [l]:", gas_consumed)
            print("Saturation of Tissue Compartments:")
            self.print_comp()
            print("****************** Descent ******************")

        # Update Dive Profile
        self.dive_profile["D"].append(str(depth_i)+" -> "+str(depth_f))
        self.dive_profile["T"].append(self.start)
        self.dive_profile["R"].append(self.run)
        self.dive_profile["S"].append("-")
        self.dive_profile["G"].append(self.gas_labels[self.current_gas])
        self.dive_profile["V"].append(gas_consumed)

        self.start = self.run

        # Now that we are at depth... Stay there
        self.constant_depth_bottom(args.tdepth, args.runT)

    def constant_depth_bottom(self, depth, time):

        # Convert depth to bar
        self.d = self.convert_to_press_abs(depth)

        # Remove time to descend to runtime
        time -= self.run

        if time <= 0:
            print("")
            print("ERROR: Insufficient runtime to reach the bottom and perform the dive.")
            exit()

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

        # Add descent to runtime
        self.run += time

        # Compute gas needed for descent
        gas_consumed = self.sac * self.d/ self.palt * time

        # Print Debug info
        if args.debug:
            print("")
            print("****************** Bottom ******************")
            print("Current depth [bar]:", self.d)
            print("Time at depth [min]:", time)
            print("")
            print("Current Runtime [min]:", self.run)
            print("Current Gas [O2/He]:", self.gas_labels[self.current_gas])
            print("Gas Volume used [l]:", gas_consumed)
            print("Saturation of Tissue Compartments:")
            self.print_comp()
            print("****************** Bottom ******************")

        # Update Dive Profile
        self.dive_profile["D"].append(str(depth))
        self.dive_profile["T"].append(self.start)
        self.dive_profile["R"].append(self.run)
        self.dive_profile["S"].append(int(self.run - self.start))
        self.dive_profile["G"].append(self.gas_labels[self.current_gas])
        self.dive_profile["V"].append(gas_consumed)


        # Once we know the saturation of the tissues, check the ascent ceiling
        self.start = self.run
        self.start_depth = self.d
        self.check_ascent_ceiling_first()

    def check_ascent_ceiling_first(self):
        """ Computes the ceiling depth given the current inert gas loading. """

        # For each compartment
        for comp_idx in range(16):

            # Compute Buhlmann M-value line parameters A and B with both N2 and He
            A = ((self.a_n2[comp_idx] * self.compartments[comp_idx,0]) + (self.a_he[comp_idx] * self.compartments[comp_idx,1]))/(self.compartments[comp_idx,2])
            B = ((self.b_n2[comp_idx] * self.compartments[comp_idx,0]) + (self.b_he[comp_idx] * self.compartments[comp_idx,1]))/(self.compartments[comp_idx,2])

            # The first time that we compute the ceiling depth: use GF_low.
            self.ascent_ceil_bar[comp_idx] = ((self.compartments[comp_idx,0]+self.compartments[comp_idx,1]) - args.gf_low * A) / ( args.gf_low/B - args.gf_low + 1)

        # Check among all tissues, which one has the highest loading.
        self.max_loading = np.max(self.ascent_ceil_bar)
                
        # Convert the highest loading to a depth increment of 3 m.
        self.current_deco_stop = self.convert_to_depth(self.max_loading)
        self.current_deco_stop = int(3 * np.ceil(self.current_deco_stop / 3))
                
        if self.current_deco_stop < 0:
            self.current_deco_stop = 0

        # Set the value of first stop given GF_low
        # This is my choice... as there are many:
        #  1) The first actual stop may not be the one computed here if we consider off-gassing in ascent.
        #  2) We could consider the increment of 3 m and off-gassing during ascent.
        #  3) We could consider the increment of 3 m and not off-gassing during ascent.
        #  4) Just self.max_loading

        self.first_stop = self.convert_to_press_abs(self.current_deco_stop)

        # Print Debug info
        if args.debug:
            print("")
            print("****************** Ceiling check FIRST ******************")
            print("Current Deco Stop [3m]:", self.current_deco_stop)
            print("Current True Stop [m]:", self.convert_to_depth(self.max_loading))
            print("Current Ascent Ceilings [bar]:", self.ascent_ceil_bar)
            print("****************** Ceiling check FIRST ******************")


        # Kill the program if there is no ceiling for the dive.
        if self.current_deco_stop == 0.0:
            print("")
            print("This is a dive within the NDL. Safety stops at:")
            print("9 m --- 1 min")
            print("6 m --- 2 min")
            print("3 m --- 3 min")
            exit()

        # Check if we need to switch gas on the way to the first ceiling.        
        gas_count = 0 
        tmp_idx = self.current_gas
        for gas_mod in self.mod_deco:
            # If the deco MOD of the gas is bigger than the first ceiling... we need to switch!
            if gas_mod - self.first_stop > 0:
                tmp_idx = gas_count                    
            gas_count += 1

        # If tmp_idx != self.current_gas... we need to switch gas on the way up.
        if tmp_idx != self.current_gas:

            swich_depth = round(self.convert_to_depth(self.mod_deco[tmp_idx])/3)*3
            self.gas_switch(swich_depth)

        # Otherwise just go to the first ceiling
        else:
            self.ascent_deep(self.current_deco_stop)

    def ascent_deep(self, depth_f):

        d1 = self.d
        d2 = self.convert_to_press_abs(depth_f)

        speed = self.speed_deep

        for comp_idx in range(16):
            p0_n2 = self.compartments[comp_idx,0]
            p0_he = self.compartments[comp_idx,1]

            pi0_n2 = (d1 - Constants.AlveolarWVP) * self.fn2[self.current_gas]
            pi0_he = (d1 - Constants.AlveolarWVP) * self.fhe[self.current_gas]

            RN2a = speed * self.fn2[self.current_gas]
            RHe = speed * self.fhe[self.current_gas]
            # The np.round is my choice ... I could leave it as a fraction but I prefer to approximate to the next integer.
            t = np.round((d2-d1) / speed)
                    
            kN2 = np.log(2) / self.ht_n2[comp_idx]
            kHe = np.log(2) / self.ht_he[comp_idx]

            pN2 = pi0_n2 + RN2a*(t - (1/kN2)) - (pi0_n2 - p0_n2 - (RN2a / kN2))* np.exp(-kN2*t)
            pHe = pi0_he + RHe*(t - (1/kHe)) - (pi0_he - p0_he - (RHe / kHe))* np.exp(-kHe*t)

            self.compartments[comp_idx,0] = pN2
            self.compartments[comp_idx,1] = pHe
            self.compartments[comp_idx,2] = pN2 + pHe

        # Add descent to runtime
        self.run += t
        
        # Compute gas needed for descent
        self.gas_consumed_ascent_deep += self.sac * ((d1+d2)/2)/ self.palt * t
        
        # Print Debug info

        if args.debug:
            print("")
            print("****************** Ascend from depth ******************")
            print("Starting depth [bar]:", d1)
            print("Final depth [bar]:", d2)
            print("Ascent time [min]:", t)
            print("Ascent speed [bar/min]:", speed)
            print("")
            print("Current Runtime [min]:", self.run)
            print("Current Gas [O2/He]:", self.gas_labels[self.current_gas])
            print("Gas Volume used [l]:", self.gas_consumed_ascent_deep)
            print("Saturation of Tissue Compartments:")
            self.print_comp()
            print("****************** Ascend from depth ******************")

        # Check if d2 is still a ceiling (check off-gassing during ascent)      
        self.d = d2
        self.check_ascent_ceiling()

        # If not go there
        if self.current_deco_stop != int(np.round(self.convert_to_depth(d2))):
            # Reset self.first_stop ?
            # self.first_stop = self.convert_to_press_abs(self.current_deco_stop)
            self.ascent_deep(self.current_deco_stop)

        else:     
            # Update dive profile
            self.dive_profile["D"].append(str(round(self.convert_to_depth(self.start_depth)))+" -> "+str(round(self.convert_to_depth(d2))))
            self.dive_profile["T"].append(self.start)
            self.dive_profile["R"].append(self.run)
            self.dive_profile["S"].append("-")
            self.dive_profile["G"].append(self.gas_labels[self.current_gas])
            self.dive_profile["V"].append(self.gas_consumed_ascent_deep)

            self.start = self.run
            self.start_depth = round(self.convert_to_depth(d2))

            # and, OK, now we are at a stable ceiling... go into the deco loop!
            self.deco()

    def check_ascent_ceiling(self):
        """ Computes the ceiling depth given the current inert gas loading. """

        # Update GF at current depth
        self.compute_GF()

        # Store the surface GFs for each compartment
        GF_surf = np.zeros(16)

        # For each compartment
        for comp_idx in range(16):

            # Compute Buhlmann M-value line parameters A and B with both N2 and He
            A = ((self.a_n2[comp_idx] * self.compartments[comp_idx,0]) + (self.a_he[comp_idx] * self.compartments[comp_idx,1]))/(self.compartments[comp_idx,2])
            B = ((self.b_n2[comp_idx] * self.compartments[comp_idx,0]) + (self.b_he[comp_idx] * self.compartments[comp_idx,1]))/(self.compartments[comp_idx,2])

            M = (self.palt / B) + A
            GF_surf[comp_idx] = (self.compartments[comp_idx,2] - self.palt) / (M -self.palt) * 100

            # Compute the ceiling depth
            self.ascent_ceil_bar[comp_idx] = ((self.compartments[comp_idx,2]) - self.GF * A) / ( self.GF/B - self.GF + 1)

        # Check among all tissues, which one has the highest loading.
        self.max_loading = np.max(self.ascent_ceil_bar)
                
        # Convert the highest loading to a depth increment of 3 m.
        self.current_deco_stop = self.convert_to_depth(self.max_loading)
        self.current_deco_stop = int(3 * np.ceil(self.current_deco_stop / 3))
                
        # We want the last positive stop to be at the user defined depth!
        if self.current_deco_stop > 0 and self.current_deco_stop < args.last_deco:
            self.current_deco_stop = args.last_deco

        if self.current_deco_stop <= 0:
            self.current_deco_stop = 0
        

        # Print Debug info
        if args.debug:          
            print("")
            print("****************** Ceiling check ******************")
            print("Current Deco Stop [3m]:", self.current_deco_stop)
            print("Current True Stop [m]:", self.convert_to_depth(self.max_loading))
            print("Current Gradient Factor [%]:", self.GF*100)
            print("Current Surface Gradient Factor [%]:", np.max(GF_surf))
            print("Current Ascent Ceilings [bar]:", self.ascent_ceil_bar)
            print("****************** Ceiling check ******************")

    def deco(self):

        if args.debug:
            print("****************** DECO PART! ******************")

        stop_min = 0
        self.gas_consumed_deco = 0

        while True:
            if self.current_deco_stop == 0:
                break
            else:
                
                # Check if there is a better gas and switch to it
                for gas_count, gas_mod in enumerate(self.mod_deco):
                    if (gas_mod + self.palt - 1 ) - self.d >= 0:
                        self.current_gas = gas_count

                # Stay 1 min at ceiling
                self.constant_depth_deco()
                stop_min += 1

                # Check deco
                self.check_ascent_ceiling()

                if self.current_deco_stop != int(np.round(self.convert_to_depth(self.d))):
                    
                    # Update dive profile
                    self.dive_profile["D"].append(str(round(self.convert_to_depth(self.d))))
                    self.dive_profile["T"].append(self.start)
                    self.dive_profile["R"].append(self.run)
                    self.dive_profile["S"].append(stop_min)
                    self.dive_profile["G"].append(self.gas_labels[self.current_gas])
                    self.dive_profile["V"].append(self.gas_consumed_deco)
                    
                    self.start = self.run
                    
                    self.ascent_shallow(self.current_deco_stop)
                    
                    stop_min = 0
                    self.start = self.run
                    self.gas_consumed_deco = 0
                else:
                    continue

    def ascent_shallow(self, depth_f):

        d1 = self.d
        d2 = self.convert_to_press_abs(depth_f)

        speed = self.speed_shallow

        for comp_idx in range(16):
            p0_n2 = self.compartments[comp_idx,0]
            p0_he = self.compartments[comp_idx,1]

            pi0_n2 = (d1 - Constants.AlveolarWVP) * self.fn2[self.current_gas]
            pi0_he = (d1 - Constants.AlveolarWVP) * self.fhe[self.current_gas]

            RN2a = speed * self.fn2[self.current_gas]
            RHe = speed * self.fhe[self.current_gas]
            t = (d2-d1) / speed
                    
            kN2 = np.log(2) / self.ht_n2[comp_idx]
            kHe = np.log(2) / self.ht_he[comp_idx]

            pN2 = pi0_n2 + RN2a*(t - (1/kN2)) - (pi0_n2 - p0_n2 - (RN2a / kN2))* np.exp(-kN2*t)
            pHe = pi0_he + RHe*(t - (1/kHe)) - (pi0_he - p0_he - (RHe / kHe))* np.exp(-kHe*t)

            self.compartments[comp_idx,0] = pN2
            self.compartments[comp_idx,1] = pHe
            self.compartments[comp_idx,2] = pN2 + pHe

        # Add descent to runtime
        self.run += t
             
        # Print Debug info

        if args.debug:
            print("")
            print("****************** Ascend from depth ******************")
            print("Starting depth [bar]:", d1)
            print("Final depth [bar]:", d2)
            print("Ascent time [min]:", t)
            print("Ascent speed [bar/min]:", speed)
            print("")
            print("Current Runtime [min]:", self.run)
            print("Current Gas [O2/He]:", self.gas_labels[self.current_gas])
            print("Gas Volume used [l]:", self.gas_consumed_ascent_deep)
            print("Saturation of Tissue Compartments:")
            self.print_comp()
            print("****************** Ascend from depth ******************")

        # Check if d2 is still a ceiling (check off-gassing during ascent)      
        self.d = d2
        if self.d > self.palt:
            self.check_ascent_ceiling()

        # Compute gas consumed
        gas_consumed = self.sac * ((d1+d2)/2)/ self.palt * t

        # Update dive profile
        self.dive_profile["D"].append(str(round(self.convert_to_depth(d1)))+" -> "+str(round(self.convert_to_depth(d2))))
        self.dive_profile["T"].append(self.start)
        self.dive_profile["R"].append(self.run)
        self.dive_profile["S"].append("-")
        self.dive_profile["G"].append(self.gas_labels[self.current_gas])
        self.dive_profile["V"].append(gas_consumed)

    def constant_depth_deco(self):

        time = 1 # one minute deco at a time for recursion

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

        # Add descent to runtime
        self.run += time
        
        # Compute gas needed for deco
        self.gas_consumed_deco += self.sac * self.d / self.palt * time

########################## Main ########################

def main():

    # Load the parameters
    params = balgo.ZHL16c_1b()
    
    # Initialize all the tissues to air partial pressures and run the mai algorithm
    tissues = Compartments(params)
    tissues.initialize()

    # Print dive profile on screen
    tissues.print_profile()

if __name__ == "__main__":
    main()
