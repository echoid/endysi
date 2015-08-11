# Script for calculating miRNA and ceRNA totals for the 10 by 10 ceRNET
# Just for testing right now; will be integrated into the main classes later

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

    # count the NaNs
    fn = join(resultsDir, ensembleName + '_allWCCs.csv')
    allCorrs = np.genfromtxt(fn, names=True, delimiter=';')
    nanCount = 0
    for i in range(len(allCorrs['r'])):
        if allCorrs['r'][i] == 0.0:
            nanCount += 1

    if nanCount > 0:
        fn = join(resultsDir, 'nanCount.txt')
        with open(fn, 'w') as nanFile:
            nanFile.write('nan count is: %d\n' % nanCount)

    # set up table for data
    dt = [('model', 'i4'), ('ceRNA_tot', 'f8'), ('miRNA_tot', 'f8'),
          ('ratio', 'f8'), ('avg_ce_trans', 'f8'), ('avg_mi_trans', 'f8'),
          ('avg_ce_deg', 'f8'), ('avg_mi_deg', 'f8'), ('avg_corr', 'f8'),
          ('min_corr', 'f8'), ('max_corr', 'f8')]

    allData = np.zeros(s, dtype=dt)

    for i in range(1, s + 1):
        modelDir = join(dataDir, 'model%d' % i)
        modelName = 'ceRNET_%dx%d_%d_%d' % (m, n, k, i)

        # Load the data for the model
        pfn = join(modelDir, modelName + ".csv")
        params = np.genfromtxt(pfn, dtype='a30', delimiter=';')

        cfn = join(modelDir, modelName + '_ceRNA_equil_WCCs.csv')
        wccs = np.genfromtxt(cfn, names=True, delimiter=';')

        dfn = join(modelDir, modelName + '_equil.cdat')
        cdat = np.genfromtxt(dfn, names=True)

        sfn = join(modelDir, modelName + '_equil_steadyStates.csv')
        ss = np.genfromtxt(sfn, names=True, delimiter=';')

        # calc steady state totals
        ceRNA_tot = 0
        miRNA_tot = 0
        for j in range(1, n + 1):
            mol = 'ceRNA%d' % j
            ceRNA_tot += ss[mol]

        for j in range(11, (2 * m) + 1):
            mol = 'S%d' % j
            miRNA_tot += cdat[mol][1]

        # calc ratio
        rnaRatio = miRNA_tot / float(ceRNA_tot)

        # calc average correlation
        avgCorr = sum(wccs['r']) / len(wccs['r'])
        minCorr = min(wccs['r'])
        maxCorr = max(wccs['r'])

        # calc avg param values
        ceTrans = 0
        ceDeg = 0
        miTrans = 0
        miDeg = 0
        for pair in params:
            if pair[0].startswith('pR'):
                ceTrans += float(pair[1])

            elif pair[0].startswith('dR'):
                ceDeg += float(pair[1])

            elif pair[0].startswith('pS'):
                miTrans += float(pair[1])

            elif pair[0].startswith('dS'):
                miDeg += float(pair[1])

        ceTrans_avg = ceTrans / float(n)
        ceDeg_avg = ceDeg / float(n)
        miTrans_avg = miTrans / float(m)
        miDeg_avg = miDeg / float(m)

        # add to the array
        allData['model'][i - 1] = i
        allData['ceRNA_tot'][i - 1] = ceRNA_tot
        allData['miRNA_tot'][i - 1] = miRNA_tot
        allData['ratio'][i - 1] = rnaRatio
        allData['avg_corr'][i - 1] = avgCorr
        allData['avg_ce_trans'][i - 1] = ceTrans_avg
        allData['avg_ce_deg'][i - 1] = ceDeg_avg
        allData['avg_mi_trans'][i - 1] = miTrans_avg
        allData['avg_mi_deg'][i - 1] = miDeg_avg
        allData['min_corr'][i - 1] = minCorr
        allData['max_corr'][i - 1] = maxCorr

    # write the data to file
    head = ''
    for name in allData.dtype.names:
        head += '%s;' % name
    head = head[:-1]

    filename = join(resultsDir, ensembleName + '_SteadyState_Analysis_Data.csv')
    np.savetxt(filename, allData, delimiter=';', header=head, comments='',
               fmt='%.18e')

    # plot results
    # totals scatter
    fig = plt.figure()
    plt.scatter(allData['ceRNA_tot'], allData['miRNA_tot'], s=10)
    fig.suptitle('Total ceRNA vs Total miRNA')
    plt.xlabel('total ceRNA')
    plt.ylabel('total miRNA')
    figName = join(resultsDir, ensembleName + '_total_ceRNA_vs_miRNA.png')
    plt.savefig(figName, dpi=120, bbox_inches='tight')
    plt.close(fig)

    # ratio vs avg correlation scatter
    fig = plt.figure()
    plt.scatter(allData['avg_corr'], allData['ratio'], s=10)
    fig.suptitle('Average correlation vs (total miRNA / total ceRNA)')
    plt.xlabel('average correlation')
    plt.ylabel('total miRNA:ceRNA ratio')
    figName = join(resultsDir, ensembleName + '_rnaRatio_vs_avgCorr.png')
    plt.savefig(figName, dpi=120, bbox_inches='tight')
    plt.close(fig)

    # min corr vs max corr
    fig = plt.figure()
    plt.scatter(allData['min_corr'], allData['max_corr'], s=10)
    fig.suptitle('Minimum correlation vs Maximum correlation')
    plt.xlabel('minimum correlation (per condition)')
    plt.ylabel('maximum correlation (per condition)')
    figName = join(resultsDir, ensembleName + '_minCorr_vs_maxCorr.png')
    plt.savefig(figName, dpi=120, bbox_inches='tight')
    plt.close(fig)

    # avg ceRNA trans vs avg miRNA trans
    fig = plt.figure()
    plt.scatter(allData['avg_ce_trans'], allData['avg_mi_trans'], s=10)
    fig.suptitle('Average ceRNA transcription rate vs average miRNA ' +
                 'transcription rate')
    plt.xlabel('average ceRNA transcription rate')
    plt.ylabel('average miRNA transcription rate')
    figName = join(resultsDir, ensembleName + '_avgCeTrans_vs_avgMiTrans.png')
    plt.savefig(figName, dpi=120, bbox_inches='tight')
    plt.close(fig)

    # avg ceRNA deg vs avg miRNA deg
    fig = plt.figure()
    plt.scatter(allData['avg_ce_deg'], allData['avg_mi_deg'], s=10)
    fig.suptitle('Average ceRNA decay rate vs average miRNA decay rate')
    plt.xlabel('average ceRNA decay rate')
    plt.ylabel('average miRNA decay rate')
    maxX = max(allData['avg_ce_deg']) + 0.001
    maxY = max(allData['avg_mi_deg']) + 0.001
    minX = min(allData['avg_ce_deg']) - 0.001
    minY = min(allData['avg_mi_deg']) - 0.001
    plt.axis([minX, maxX, minY, maxY])
    figName = join(resultsDir, ensembleName + '_avgCeDeg_vs_avgMiDeg.png')
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
    parser.add_argument('--time', type=str,
                        help='timestamp for the run')
    parser.add_argument('--root', type=str, default='None',
                        help='timestamp for the run')
    args = parser.parse_args()

    if args.root == 'None':
        doit(args.m, args.n, args.k, args.s, args.time)
    else:
        doit(args.m, args.n, args.k, args.s, args.time, rootDir=args.root)