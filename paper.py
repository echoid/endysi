#!/usr/bin/env python

###############################################################################
# This script runs the experiments and generates the figures for the paper
# to be published.
#
###############################################################################

import socket
import numpy as np
import matplotlib.pyplot as plt
import endysi
import parameters as params
from utilities import *


def _writeDataToCSV(filename, data, delim=';', comment='', frmt='%.18e'):
    head = ''
    for name in data.dtype.names:
        head += '%s;' % name
    head = head[:-1]

    np.savetxt(filename, data, delimiter=delim, header=head,
               comments=comment, fmt=frmt)


def _plotData_semilogx(filename, data, xName, xLabel, yLabel, ddpi=120):
    names = [name for name in data.dtype.names if name != xName]

    fig = plt.figure()
    for name in names:
        plt.semilogx(data[xName], data[name], label=name)

    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.legend()
    plt.savefig(filename, dpi=ddpi, bbox_inches='tight')
    plt.close(fig)


def figure1():
    paramRange = params.Figure1Parameters()

    m = 1
    n = 2
    k = 1
    s = 1000
    outFreq = 10000
    method = 'ssa'
    maxHalfLife = 70000
    halfLifeMults = 2
    nSamples = 1000

    tEnd = maxHalfLife * halfLifeMults * nSamples

    # Check where we're running
    baseDir = None
    host = socket.gethostname()
    if host == 'crick':
        baseDir = '/ohri/projects/perkins/mattm/ceRNA/endysi'
    elif host == 'ogic':
        baseDir = '/data/matt/ceRNA/endysi'

    # make the ensemble and run it
    ens = endysi.Ensemble(m, n, k, s, method, tEnd, outFreq, nSamples,
                   paramRange, baseDir=baseDir, name='Figure1')
    ens.runAll()

    # Plot additional results
    # Make a table for the data
    dt = [('miRNA trans rate', 'f8'), ('norm ceRNA1', 'f8'),
          ('norm ceRNA2', 'f8'), ('norm miRNA', 'f8'), ('pearson r', 'f8')]

    data = np.zeros(s, dtype=dt)

    # Gather the data
    data['miRNA trans rate'] = ens.pRange

    # Normalize steady states
    for i in range(s):
        ce1 = ens.equilSS['ceRNA1'][i]
        ce2 = ens.equilSS['ceRNA2'][i]
        mi = ens.equilSS['miRNA1'][i]
        denom = ce1 + ce2 + mi
        normR1 = ce1 / denom
        normR2 = ce2 / denom
        normS = mi / denom
        data['norm ceRNA1'][i] = normR1
        data['norm ceRNA2'][i] = normR2
        data['norm miRNA'][i] = normS

    data['pearson r'] = ens.ceWCCCs

    # Write data to file
    fn = join(ens.resultsDir, ens.name + '_normSS_and_R_data.csv')
    _writeDataToCSV(fn, data)

    # Plot data
    fn = join(ens.resultsDir, ens.name + '_normSS_and_R.png')
    _plotData_semilogx(fn, data, 'miRNA trans rate',
                       'miRNA transcription rate', '')


def figure2(part):
    # Check where we're running
    baseDir = None
    host = socket.gethostname()
    if host == 'crick':
        baseDir = '/ohri/projects/perkins/mattm/ceRNA/endysi'
    elif host == 'ogic':
        baseDir = '/data/matt/ceRNA/endysi'

    if part == 1:
        m = 1
        n = 2
        k = 1
        s = 100000
        outFreq = 5000
        method = 'ode'
        maxHalfLife = 70000
        halfLifeMults = 2
        nSamples = 1

        tEnd = maxHalfLife * halfLifeMults * nSamples

        corrs = []
        pMin = 1.0e-6
        pMax = 1.0e-4
        resultsDir = None
        for i in range(5):
            paramRange = params.NitzanParametersCustom(pMin, pMax)
            pMin *= 10.0
            pMax *= 10.0

            # make the ensemble and run it
            ens = endysi.Ensemble(m, n, k, s, method, tEnd, outFreq, nSamples,
                                  paramRange, baseDir=baseDir, name='Figure2.1')
            ens.runAll()

            corrs.append(np.mean(ens.ceEquilCCs['r']))
            if resultsDir is None:
                resultsDir = ens.rRunDir

        fn = join(resultsDir, 'Figure2_ceCCCs_%s.csv' % genTimeString())
        np.savetxt(fn, np.array(corrs))

    if part == 2:
        m = 1
        n = 2
        k = 1
        s = 1000
        outFreq = 10000
        method = 'ssa'
        maxHalfLife = 70000
        halfLifeMults = 2
        nSamples = 1000

        tEnd = maxHalfLife * halfLifeMults * nSamples

        corrs = []
        pMin = 1.0e-6
        pMax = 1.0e-4
        resultsDir = None
        for i in range(5):
            paramRange = params.NitzanParametersCustom(pMin, pMax)
            pMin *= 10.0
            pMax *= 10.0

            # make the ensemble and run it
            ens = endysi.Ensemble(m, n, k, s, method, tEnd, outFreq, nSamples,
                                  paramRange, baseDir=baseDir, name='Figure2.2')
            ens.runAll()

            corrs.append(np.mean(ens.ceWCCCs))
            if resultsDir is None:
                resultsDir = ens.rRunDir

        fn = join(resultsDir, 'Figure2_ceWCCs_%s.csv' % genTimeString())
        np.savetxt(fn, np.array(corrs))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--figure', type=str, help='The figure to generate')

    args = parser.parse_args()

    if args.figure == '1':
        figure1()
    elif args.figure == '2p1':
        figure2(1)
    elif args.figure == '2p2':
        figure2(2)
