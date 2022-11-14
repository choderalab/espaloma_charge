import sys
import argparse
from espaloma_charge.antechamber_utils import charge_mol2


def espaloma_charge():
    parser = argparse.ArgumentParser(description="Espaloma charge.")
    parser.add_argument("-i", type=str, required=True)
    parser.add_argument("-o", type=str, default="")
    args = parser.parse_args()
    if ".mol2" not in args.i:
        raise NotImplementedError
    charge_mol2(args.i, args.o)


if __name__ == "__main__":
    espaloma_charge()
