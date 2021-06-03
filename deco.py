import numpy as np
import matplotlib.pyplot as plt

# Important note: all pressure units are expressed in [Bar].

############# USER-DEFINED VARIABLES ################
repetitive = False
h_on_sea = 0.0
last_deco_stop = 3 # in meters

fN2=0.7808
fHe=0.0


pW=0.0626615132
atm_bar = 1.01325

################### CONSTANTS #######################

WV_PRESSURE = 0.0627            # water vapor pressure in bar, based on respiratory quotient Rq = 1.0 (Buhlmann value)
WV_PRESSURE_SCHREINER = 0.0493  # water vapor pressure in bar, based on respiratory quotient Rq = 0.8 (Schreiner value)

#####################################################

def corr_altitude_alveoli(fN2=0.7808, pa=1.01325, pW=0.0626615132):
    """ Corrects nitrogen partial pressure for altitude and alveolar water vapour.
    fN2: Fraction of N2 in breathed gas
    pa = Atmospheric pressure. (sea level = 1.01325 bar)
    pw = Water vapour pressure in alveoli """
    
    ppN2 = fN2 * ( pa - pW )

    return ppN2

def altitude(h_on_sea=0.0):
    """ Returns the total pressure of air given the height over sea-level in meters."""
    
    # Constants for the barometric formula

    _M = 0.02896 # Molar mass of earth's air
    _g = 9.807 # Earth's gravitational acceleration
    _R = 8.3143 # Ideal gas constant
    _T = 288.15 # Standard temperature
    _P0 = 1.01325 # Pressure at sea-level in [Bar]

    cnst = (_M*_g)/(_R*_T)

    return _P0 * np.exp( -cnst * h_on_sea )

def initialize():
    """ Initialize inert gasses partial pressures and total inert gas pressure for all the 16 tissues."""
    palt = altitude(h_on_sea)

    pN2 = corr_altitude_alveoli(pa=palt)
    pHe = 0
    pTot = pN2 + pHe
    
    tissues = np.zeros((16,3))
    for i in range(16):
        tissues[i,0] = pN2
        tissues[i,1] = pHe
        tissues[i,2] = pTot

    return tissues

def constant_ascent_descent(idepth, fdepth, time, tissues, ht):
    """Returns the tissue loading when the diver ascends or descends at constant speed.
    Please not that descent speed is positive and ascent is negative."""

    # Convert depth from meters to pressure in [Bar]
    idepth = ((idepth / 10) + 1) * atm_bar
    fdepth = ((fdepth / 10) + 1) * atm_bar

    speed = (fdepth-idepth) / time

    for i in range(16):
        
        P0_N2 = tissues[i,0]
        P0_He = tissues[i,1]

        Pdepth_N2 = (idepth - pW) * fN2
        Pdepth_He = (idepth - pW) * fHe

        RN2 = speed * fN2
        RHe = speed * fHe
        
        Kn = np.log(2) / ht[i,0]
        Khe = np.log(2) / ht[i,1]

        Pn = RN2/Kn*(np.exp(-Kn*time) + Kn*time -1) +  P0_N2 * np.exp(-Kn*time) + Pdepth_N2 * (1 - np.exp(-Kn*time))
        Phe = RHe/Khe*(np.exp(-Khe*time) + Khe*time -1) +  P0_He * np.exp(-Khe*time) + Pdepth_He * (1 - np.exp(-Khe*time))

        tissues[i,0] = Pn
        tissues[i,1] = Phe

        tissues[i,2] = Pn + Phe

    return tissues

def constant_depth(depth, t, tissues, ht):
    """ Returns the tissue loading when the diver reamins at constant depth."""

    # Convert depth from meters to pressure in [Bar]
    depth = ((depth / 10) + 1) * atm_bar

    for i in range(16):

        P0_N2 = tissues[i,0]
        P0_He = tissues[i,1]

        Pdepth_N2 = (depth - pW) * fN2
        Pdepth_He = (depth - pW) * fHe

        Pn = P0_N2 * 2**(-t / ht[i,0]) + Pdepth_N2 * (1 - 2**(-t / ht[i,0]))
        Phe = P0_He * 2**(-t / ht[i,1]) + Pdepth_He * (1 - 2**(-t / ht[i,1]))

        tissues[i,0] = Pn
        tissues[i,1] = Phe

        tissues[i,2] = Pn + Phe

    return tissues

def check_ceiling(tissues, ht):
    """ Compute the ceiling depth based on the Buhlmann A and B values and the inert gasses partial pressure """

    all_ceil = []
    for i in range(16):
    # Compute the A and B value as a weighted average of the A and B value of each inert gas
        A = ((ht[i,2] * tissues[i,0]) + (ht[i,4] * tissues[i,1]))/tissues[i,2]
        B = ((ht[i,3] * tissues[i,0]) + (ht[i,5] * tissues[i,1]))/tissues[i,2]
        ceil = (tissues[i,2] - A) * B
        ceil_in_m = (ceil - atm_bar)*10

        # Round to next multiple of 3m
        all_ceil.append(3 * round(ceil_in_m/3))

    all_ceil = np.array(all_ceil)

    if max(all_ceil) <= 0:
        return 0
    else:
        return max(all_ceil)

def length_of_deco(curr_stop, tissues, ht):

    # Compute all stops in steps of 3m
    all_stops = np.arange(curr_stop, last_deco_stop -3, -3)

    # Convert stop length in Bar
    all_stops = ((all_stops / 10) + 1) * atm_bar
    all_stops = np.append(all_stops, atm_bar)

    for i,stop in enumerate(all_stops[:-1]):
        stop_length = 0 
        while((check_ceiling(tissues, ht) / 10 + 1) * atm_bar > all_stops[i+1]):
            # Stay 1 min and update the inert gas partial pressure
            tissues = constant_depth(stop, 1.0, tissues, ht)
            stop_length += 1

        print(np.floor((stop/atm_bar -1 )*10) , stop_length)

    return 0

def no_deco():
    return 0

def main():

    # Import half_times
    ht = np.loadtxt("half_times_ZH-L16c.csv", comments='#', delimiter=",")
    
    # If the diver spent enough time on surface.
    if not repetitive:
        # Initialize all the tissues to air partial pressures.
        tissues = initialize()

    tissues = constant_ascent_descent(0, 40, 15, tissues, ht)
    tissues = constant_depth(40, 20, tissues, ht) 
    ceil = check_ceiling(tissues, ht)

    if ceil == 0:
        print("No Deco")
    else:
        print("Ceiling: "+str(ceil)+" m")
        length_of_deco(ceil, tissues, ht)

    tissues = constant_ascent_descent(40, 15, 15, tissues, ht)
    ceil = check_ceiling(tissues, ht)

    if ceil == 0:
        print("No Deco")
    else:
        print("Ceiling: "+str(ceil)+" m")
        length_of_deco(ceil, tissues, ht)






if __name__ == "__main__":
    main()