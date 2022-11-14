import sys
import argparse

def espaloma_charge():
    parser = argparse.ArgumentParser(description="Espaloma charge.")
    parser.add_argument("-i", type=str)
    parser.add_argument("-o", type=str, default="")
    args = parser.parse_args()
    if ".mol2" not in args.i:
        raise NotImplementedError
    from .antechamber_utils import charge_mol2
    charge_mol2(args.i, args.o)
