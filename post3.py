# Script for steady state analysis of 10 by 10 ceRNA networks.
# Just for testing right now; will be integrated into the main classes later.
# ODEs, CCCCs, 100 sets, 100 repeats

from __future__ import print_function, division
import os
import numpy as np
import matplotlib.pyplot as plt
from utilities import *

# global attributes
nModels = 100


def doit(m, n, k, s, p, timestamp, rootDir=None):
    # Set up data directories
    popName = 'Pop@%d_%dx%dc%dx%d' % (p, m, n, k, s)
    if rootDir is None:
        rootDir = join(os.path.expanduser('~'), 'research/results/ceRNA/endysi')

    homeDir = join(rootDir, popName + '/' + timestamp)
    dataDir = join(homeDir, 'data')
    ensName = 'ceRNET_%dx%dc%dx%d' % (m, n, k, s)
    netDir = join(dataDir, ensName)
    resultsDir = join(homeDir, 'results')

    # set up table for data
    dt = [('ens', 'i4'), ('avg_ceRNA', 'f8'), ('avg_miRNA', 'f8'),
          ('ratio', 'f8'), ('avg_corr', 'f8')]

    allData = np.zeros(p, dtype=dt)

    for i in range(1, p + 1):
        ensDir = join(netDir, 'e%d' % i)
        ensResultsDir = join(ensDir, 'results')
        ssTabName = join(ensResultsDir, ensName + '_steadyStates.csv')
        ssTab = np.genfromtxt(ssTabName, delimiter=';', names=True)
        ceRNAs = [name for name in ssTab.dtype.names if name.startswith('ce')]
        miRNAs = [name for name in ssTab.dtype.names if name.startswith('mi')]
        ceMeans = []
        miMeans = []
        for ceRNA in ceRNAs:
            ceMeans.append(np.mean(ssTab[ceRNA]))
        for miRNA in miRNAs:
            miMeans.append(np.mean(ssTab[miRNA]))

        ceMom = np.mean(ceMeans)
        miMom = np.mean(miMeans)

        ccTabName = join(ensResultsDir, ensName + '_ceRNA_CCs.csv')
        ccTab = np.genfromtxt(ccTabName, names=True, delimiter=';')
        avgR = np.mean(ccTab['r'])

        allData['ens'][i - 1] = i
        allData['avg_ceRNA'][i - 1] = ceMom
        allData['avg_miRNA'][i - 1] = miMom
        allData['avg_corr'][i - 1] = avgR
        allData['ratio'][i - 1] = ceMom / miMom

    # write the data to file
    head = ''
    for name in allData.dtype.names:
        head += '%s;' % name
    head = head[:-1]

    filename = join(resultsDir, popName + '_SteadyState_Analysis_Data.csv')
    np.savetxt(filename, allData, delimiter=';', header=head, comments='',
               fmt='%.18e')

    # plot results
    # totals scatter
    fig = plt.figure()
    plt.scatter(allData['avg_ceRNA'], allData['avg_miRNA'], s=10)
    fig.suptitle('average ceRNA vs average miRNA')
    plt.xlabel('Avg ceRNA')
    plt.ylabel('Avg miRNA')
    figName = join(resultsDir, popName + '_avg_ceRNA_vs_miRNA.png')
    plt.savefig(figName, dpi=120, bbox_inches='tight')
    plt.close(fig)

    # ratio vs avg correlation scatter
    fig = plt.figure()
    plt.scatter(allData['avg_corr'], allData['ratio'], s=10)
    fig.suptitle('Average correlation vs (avg miRNA / total ceRNA)')
    plt.xlabel('average correlation')
    plt.ylabel('avg miRNA:ceRNA ratio')
    figName = join(resultsDir, popName + '_rnaRatio_vs_avgCorr.png')
    plt.savefig(figName, dpi=120, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=int, help='The number of miRNAs')
    parser.add_argument('-n', type=int, help='The number of ceRNAs')
    parser.add_argument('-k', type=int, help='The degree of each ceRNA')
    parser.add_argument('-s', type=int,
                        help='The number of models in each ensemble')
    parser.add_argument('-p', type=int,
                        help='The number of ensembles in each population')
    parser.add_argument('--time', type=str,
                        help='timestamp for the run')
    parser.add_argument('--root', type=str, default='None',
                        help='timestamp for the run')
    args = parser.parse_args()

    if args.root == 'None':
        doit(args.m, args.n, args.k, args.s, args.p, args.time)
    elif args.root == 'ohri':
        r = '/ohri/projects/perkins/mattm/ceRNA/endysi'
        doit(args.m, args.n, args.k, args.s, args.p, args.time, rootDir=r)