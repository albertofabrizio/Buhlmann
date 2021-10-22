import argparse
import numpy as np

########################## Parsing user defined input ##########################
parser = argparse.ArgumentParser(description='This simple program computes the final mix, M.O.D, E.A.D., and the E.N.D. given a starting blend and the additional gasses.')

parser.add_argument('--init_press', type=int, dest='p_init', default=0, help="Initial pressure in the tank [bar].")
parser.add_argument('--ifo2', type=float, dest='fo2', default=0.0, help="Initial fraction of oxygen in the tank.")
parser.add_argument('--ifhe', type=float, dest='fhe', default=0.0, help="Initial fraction of helium in the tank.")

parser.add_argument('--addfo2', type=float, dest='afo2', required=True, nargs='+', help="Fraction(s) of oxygen in the added gas(es).")
parser.add_argument('--addfhe', type=float, dest='afhe', required=True, nargs='+', help="Fraction(s) of helium in the added gas(es).")
parser.add_argument('--add_press', type=int, dest='p_add', required=True, nargs='+', help="Pressure(s) of the added gas(es)[bar].")

parser.add_argument('--mod', type=float, dest='mod', default=1.4, help="User defined ppO2 at maximum operating depth [bar].")

parser.add_argument('--tdepth', type=int, dest='tdepth', required=True, help="Target depth below surface [m].")



args = parser.parse_args()

########################## Helper Functions ########################

def compute_final_blend(pi, ifo2, ifhe, pa, afo2, afhe):
    """ Takes gas inputs and return the final composition of the blend."""
    # Note the compressibility of the gasses, the thermal expansion and other real gas behavior is not taken into account.
    
    # Convert to array
    pa = np.array(pa)
    afo2 = np.array(afo2)
    afhe = np.array(afhe)

    # Compute the partial pressures of each gas
    ffo2 = (ifo2 * pi + np.dot(pa, afo2)) / (pi + pa.sum())
    ffhe = (ifhe * pi + np.dot(pa, afhe)) / (pi + pa.sum())
    ffn2 = 1 - ffo2 - ffhe

    # Print on screen and return
    print("Gas fractions in the final blend.")
    print("O2:", np.round(ffo2*100,2),"%")
    print("He:", np.round(ffhe*100,2),"%")
    print("N2:", np.round(ffn2*100,2),"%")

    return ffo2, ffhe, ffn2

def compute_MOD(blend):
    """ Computes the maximum operating depth (MOD) given the blend."""

    mod = ((args.mod / blend[0]) - 1) * 10
    print("Maximum Operating Depth [m]:", np.round(mod,1))
    
    return mod

def compute_END(blend, tod):
    """ Computes the equivalent narcotic depth (END) given the blend."""

    end = (((1 - blend[1]) * tod) - 1) * 10
    print("Equivalent Narcotic Depth [m]:", np.round(end,1))
    
    return end

def compute_EAD(blend, tod):
    """ Computes the equivalent air depth (EAD) given the blend."""
    
    air_nitrogen = 0.79 # Approximate, but consistent for computation
    ead = (((blend[2] / air_nitrogen) * tod) - 1) * 10
    print("Equivalent Air Depth [m]:", np.round(ead,1))
    
    return ead

########################## Main ########################

def main():
    
    if args.fo2 == 0.0 and args.fhe == 0.0 and args.p_init !=0:
        print("Warning: Initial tank is filled up to ",args.p_init,"bar, with pure N2.")

    # Compute the final blend
    blend = np.zeros(3)
    blend = compute_final_blend(args.p_init, args.fo2, args.fhe, args.p_add, args.afo2, args.afhe)
    print("")

    # Compute MOD
    compute_MOD(blend)

    # Convert depth to Absolute Pressure
    tod = (args.tdepth / 10) + 1 # The +1 approximation is consistent with the functions END and EAD

    # Compute EAD
    compute_EAD(blend, tod)

    # Compute END (Oxygen is narcotic)
    compute_END(blend, tod)
    print("")

if __name__ == "__main__":
    main()
