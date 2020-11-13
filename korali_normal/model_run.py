#!/usr/bin/env python3

import sys
import argparse

sys.path.append('./_model')
from model_call import RESULT_TAG
from model import model
from model_data import getReferencePoints

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('beta', type=float)
    parser.add_argument('sigma', type=float)
    args = parser.parse_args()

    s = dict()
    s['Parameters'] = [args.beta, args.sigma]
    model(s, getReferencePoints())

    print(RESULT_TAG)
    keys = ["Reference Evaluations", "Standard Deviation"]
    for k in keys:
        print("'{}' : {},".format(k, s[k]))

if __name__ == "__main__":
    main()
