# Another script for steady state analysis of 10 by 10 ceRNA networks.
# Here, we're interested in how many of the networks have all ceRNAs with
# steady states below 1 (and conversely, how many have all miRNAs above 1).
# For testing and validation; will be integrated into the main classes later.

from __future__ import print_function, division
import os
import numpy as np
import matplotlib.pyplot as plt
from utilities import *

# global attributes
nModels = 100


def doit(m, n, k, s, timestamp, rootDir=None):
    # Set up data directories
    ensembleName = 'ceRNET_%dx%dc%dx%d' % (m, n, k, s)
    if rootDir is None:
        rootDir = join(os.path.expanduser('~'), 'research/results/ceRNA/endysi')

    homeDir = join(rootDir, ensembleName + '/' + timestamp)
    dataDir = join(homeDir, 'data')
    resultsDir = join(homeDir, 'results')

    # Load steady state data
    fn = join(resultsDir, ensembleName + '_equil_steadyStates.csv')
    ssData = np.genfromtxt(fn, names=True, delimiter=';')

    # Add up number of networks with ceRNAs below 1 and miRNAs above 1
    ceRNAcount = 0
    miRNAcount = 0
    ceRNAs = [name for name in ssData.dtype.names if 'ceRNA' in name]
    miRNAs = [name for name in ssData.dtype.names if 'miRNA' in name]
    goodNetworks = []

    for i in range(s):
        # ceRNAs
        count = 0
        for ceRNA in ceRNAs:
            if ssData[ceRNA][i] < 1.0:
                count += 1

        if count == 10:
            ceRNAcount += 1
        else:
            goodNetworks.append(i + 1)

        # miRNAs
        count = 0
        for miRNA in miRNAs:
            if ssData[miRNA][i] > 1.0:
                count += 1

        if count == 10:
            miRNAcount += 1

    print('Number of networks with ceRNAs below 1 is', ceRNAcount)
    print('Number of networks with miRNAs above 1 is', miRNAcount)
    if len(goodNetworks) > 0:
        print(goodNetworks)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=int, help='The number of miRNAs')
    parser.add_argument('-n', type=int, help='The number of ceRNAs')
    parser.add_argument('-k', type=int, help='The degree of each ceRNA')
    parser.add_argument('-s', type=int,
                        help='The number of models in each ensemble')
    parser.add_argument('--time', type=str,
                        help='timestamp for the run')
    parser.add_argument('--root', type=str, default='None',
                        help='root directory for the run')
    args = parser.parse_args()

    if args.root == 'None':
        doit(args.m, args.n, args.k, args.s, args.time)
    else:
        doit(args.m, args.n, args.k, args.s, args.time, rootDir=args.root)