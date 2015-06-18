#!/usr/bin/env python2

from __future__ import print_function, division
import os
import math
import time
import random
import logging
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import bngl
from utilities import *

# Global debugging level
_debugLevel = logging.INFO  # use logging.INFO to turn off debugging

# Global logging switch
_logging = False

# Global dpi for plots
_dpi = 120


# Global functions (mostly for analysis and plotting)
def _plotAllTrajectories(filename, data, mols, colours=None, ddpi=120):
    """Plot the trajectories of a set of molecules.  Assumes that the
    names in the mols argument are of the form 'mol_free', as in the
    output files from the BNG simulation.  Also assumes that the
    filename argument is the full absolute path for the file."""

    if len(mols) == 0:
        return

    if colours is None:
        # colours
        purple = '#673594'
        magenta = '#a0228d'
        pink = '#ed1991'
        red = '#ee393b'
        yellow = '#f4ee37'
        green = '#4bbc44'
        blue = '#30c6f5'
        navy = '#33689c'
        aqua = '#30f5bb'
        orange = '#f58d30'

        colours = [purple, magenta, pink, red, yellow, green, blue,
                   navy, aqua, orange]

    count = 0
    fig = plt.figure()
    for mol in mols:
        plt.plot(data['time'], data[mol], colours[count],
                 label=mol.split('_')[0])
        count += 1

    plt.xlabel('Time (s)')
    plt.ylabel('Free molecules')
    plt.legend()
    plt.savefig(filename, dpi=ddpi, bbox_inches='tight')
    plt.close(fig)


def _histogram(filepath, x, xLabel, nBins=20, ddpi=120):
    if len(x) == 0:
        return

    # Check to make sure filepath doesn't include an extension
    if '.' in filepath:
        i = filepath.index('.')
        filepath = filepath[:i]

    # Make the figure
    fig = plt.figure()
    (counts, bins, patches) = plt.hist(x, bins=nBins)
    plt.xlabel(xLabel)
    plt.ylabel('frequency')
    plt.savefig(filepath + '.png', dpi=ddpi, bbox_inches='tight')
    plt.close(fig)

    # Save the data
    histData = np.zeros(nBins, dtype=[('bin', 'f8'), ('count', 'f8')])
    np.copyto(histData['bin'], bins[:len(counts)])
    np.copyto(histData['count'], counts)
    _writeDataToCSV(filepath + '.csv', histData)


def _plotAndSave(filename, x, y, xLabel, yLabel, ddpi=120):
    pass


def _writeDataToCSV(filename, data, delim=';', comment='', frmt='%.18e'):
    head = ''
    for name in data.dtype.names:
        head += '%s;' % name
    head = head[:-1]

    np.savetxt(filename, data, delimiter=delim, header=head,
               comments=comment, fmt=frmt)


### Class definitions ###
class Experiment:
    def __init__(self, method, model, tEnd, outFreq, pTarget='ceRNA', pProp=5.0):

        self.model = model
        self.simulator = bngl.BnglSimulator(model)
        self.perturbTarget = pTarget
        self.perturbProp = pProp

        self.action = 'simulate({method=>"%(meth)s",suffix=>"%(suf)s",' + \
                      'continue=>%(cnt)d,steady_state=>%(ss)d,' + \
                      't_end=>%(tend)d,n_steps=>%(nsteps)d,' + \
                      'print_CDAT=>%(pcdat)d})'

        nsteps = int(tEnd / outFreq)
        self.opts = {'suf': 'sim', 'tend': tEnd, 'nsteps': nsteps, 'cnt': 0,
                     'meth': method, 'pcdat': 0, 'ss': 0}
        return

    def run(self):
        # Initialize simulator
        os.chdir(self.model.home)
        self.simulator.initialize()

        # Equilibration phase
        self.opts['suf'] = 'equil'
        self.simulator.sendAction(self.action % self.opts)
        self.simulator.saveConcentrations()

        # Perturbation phase
        #self.opts['suf'] = 'perturb'
        #pMol = None
        #param = None

        #if self.perturbTarget == 'miRNA':
            #pMol = self.model.getMolWithName('miRNA1')
            #param = 'pR_1'
        #else:
            #pMol = self.model.getMolWithName('ceRNA1')
            #param = 'pT_1'

        #oldProdRate = pMol.prodRate
        #newProdRate = oldProdRate * self.perturbProp
        #self.simulator.setParameter(param, newProdRate)
        #self.simulator.sendAction(self.action % self.opts)
        self.simulator.close()
        del self.simulator
        return

    def loadData(self):
        # Load data from the two simulations
        self.equilData = np.genfromtxt('%s_equil.gdat' % self.model.filePath,
                                       names=True)
        #self.perturbData = np.genfromtxt('%s_perturb.gdat' %
        #                                 self.model.filePath, names=True)

        # Get and save molecule names
        self.mols = [name for name in self.equilData.dtype.names if
                     name != 'time']
        self.miRNAs = [mol for mol in self.mols if 'miRNA' in mol]
        self.ceRNAs = [mol for mol in self.mols if 'ceRNA' in mol]

        # Pair up molecules for correlations
        self.corrMols = {}
        self.corrMols['ceRNA'] = list(combinations(self.ceRNAs, 2))
        self.corrMols['miRNA'] = list(combinations(self.miRNAs, 2))
        return

    def calcSteadyStates(self):
        # Create a table for the steady state values
        ssdt = [(mol.split('_')[0], 'f8') for mol in self.mols]
        self.equilSS = np.zeros(1, dtype=ssdt)
        #self.perturbSS = np.zeros(1, dtype=ssdt)

        # Calculate steady states
        for mol in self.mols:
            molName = mol.split('_')[0]
            if self.opts['meth'] == 'ssa':
                sp = len(self.equilData['time']) - 1000
                self.equilSS[molName] = np.mean(self.equilData[mol][sp:])
                #self.perturbSS[molName] = np.mean(self.perturbData[mol][sp:])
            else:
                self.equilSS[molName] = self.equilData[mol][-1]
                #self.perturbSS[molName] = self.perturbData[mol][-1]

        # Write steady states to file
        fn = self.model.filePath + '_equil_steadyStates.csv'
        _writeDataToCSV(fn, self.equilSS)

        #fn = self.model.filePath + '_perturb_steadyStates.csv'
        #_writeDataToCSV(fn, self.perturbSS)
        return

    def plotTrajectories(self, ddpi=120):
        # Plot miRNAs
        fn = join(self.model.plotDir, self.model.name + '_miRNA_equil_traj')
        _plotAllTrajectories(fn, self.equilData, self.miRNAs)

        #fn = join(self.model.plotDir, self.model.name + '_miRNA_perturb_traj')
        #_plotAllTrajectories(fn, self.perturbData, self.miRNAs)

        # Plot ceRNAs
        fn = join(self.model.plotDir, self.model.name + '_ceRNA_equil_traj')
        _plotAllTrajectories(fn, self.equilData, self.ceRNAs)

        #fn = join(self.model.plotDir, self.model.name + '_ceRNA_perturb_traj')
        #_plotAllTrajectories(fn, self.perturbData, self.ceRNAs)
        return

    def makeFoldPlots(self, ddpi=120):
        # Create colours for the plots
        colours = []
        p = [str(x) for x in range(10)]
        p.extend(['a', 'b', 'c', 'd', 'e', 'f'])

        while len(colours) < len(self.corrMols['ceRNA']):
            s = random.sample(p, 6)
            colour = '#' + ''.join(s)
            if colour != '#ffffff' and colour not in colours:
                colours.append(colour)

        # Equil data
        i = 0
        fig = plt.figure()
        for pair in self.corrMols['ceRNA']:
            plt.plot(self.equilData[pair[0]], self.equilData[pair[1]],
                     colours[i])
            i += 1
        fn = join(self.model.plotDir, self.model.name + '_equil_foldPlots.png')
        plt.savefig(fn, dpi=ddpi, bbox_inches='tight')
        plt.close(fig)

        # Perturb data
        #i = 0
        #fig = plt.figure()
        #for pair in self.corrMols['ceRNA']:
            #plt.plot(self.perturbData[pair[0]], self.perturbData[pair[1]],
                     #colours[i])
            #i += 1
        #fn = join(self.model.plotDir, self.model.name + '_perturb_foldPlots.png')
        #plt.savefig(fn, dpi=ddpi, bbox_inches='tight')
        #plt.close(fig)
        return

    def runAnalyses(self):
        self.loadData()
        self.calcSteadyStates()
        self.plotTrajectories()
        #self.makeFoldPlots()
        return

    def deleteData(self):
        del self.equilData
        #del self.perturbData
        return


class Ensemble:
    def __init__(self, m, n, k, size, method, tEnd, outFreq, paramDict,
                 timestamp=None, baseDir=None):

        if timestamp is None:
            self.timestamp = genTimeString()
        else:
            self.timestamp = timestamp

        self.name = 'ceRNET_%dx%dc%dx%d' % (m, n, k, size)
        self.m = m
        self.n = n
        self.k = k
        self.size = size
        self.models = []
        self.experiments = []

        self.createDirectories(baseDir)

        if _logging:
            self.createLoggers()
            msg = 'Initializing system to run {0} networks with M={1}, N={2}, K={3}'
            self.logger.info(msg.format(size, m, n, k))
            self.logger.info('Parameter ranges are as follows:')
            self.logger.info(str(paramDict))

            self.logger.info('Creating models...')

        self.createModels(paramDict)
        self.createExperiments(method, tEnd, outFreq)
        return

    def createLoggers(self):
        # create loggers and set levels
        self.debugger = logging.getLogger('debugger')
        self.debugger.setLevel(_debugLevel)
        self.logger = logging.getLogger('logger')
        self.logger.setLevel(logging.INFO)

        # create file handler which writes important events to a log file
        filename = join(self.curRun, 'run_log.log')
        fh = logging.FileHandler(filename, mode='w')
        fh.setLevel(logging.INFO)

        # create console handler for writing debugging messages to the console
        ch = logging.StreamHandler()
        ch.setLevel(_debugLevel)

        # create formatter and add it to the handlers
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        # add the handlers to the logger
        self.logger.addHandler(fh)
        self.debugger.addHandler(ch)
        return

    def createDirectories(self, baseDir):
        if baseDir is None:
            self.rootDir = join(os.path.expanduser('~'),
                                'research/results/ceRNA/endysi/' + self.name)
        else:
            self.rootDir = join(baseDir, self.name)

        self.curRun = join(self.rootDir, self.timestamp)
        self.dataDir = join(self.curRun, 'data')
        self.resultsDir = join(self.curRun, 'results')
        makeDirs(self.dataDir)
        makeDirs(self.resultsDir)
        return

    def createModels(self, paramDict):
        for i in range(1, self.size + 1):
            dDir = join(self.dataDir, 'model%d' % i)
            model = bngl.CernetModel(dDir, self.m, self.n, self.k, i,
                                     paramDict, seed=None)

            self.models.append(model)
        return

    def createExperiments(self, method, tEnd, outFreq):
        for model in self.models:
            e = Experiment(method, model, tEnd, outFreq)
            self.experiments.append(e)

    def runExperiments(self):
        for experiment in self.experiments:
            if _logging:
                self.logger.info('Running simulations on %s' % experiment.model.name)

            experiment.run()
            experiment.runAnalyses()
            experiment.deleteData()
        return

    def collectSteadyStates(self):
        # Create a table for the steady state values
        names = self.experiments[0].equilSS.dtype.names
        ssdt = [(name, 'f8') for name in names]
        ssdt.insert(0, ('model', 'i4'))
        self.equilSS = np.zeros(self.size, dtype=ssdt)
        #self.perturbSS = np.zeros(self.size, dtype=ssdt)

        for experiment in self.experiments:
            i = experiment.model.index
            self.equilSS['model'][i - 1] = i
            #self.perturbSS['model'][i - 1] = i
            for name in names:
                self.equilSS[name][i - 1] = experiment.equilSS[name][0]
                #self.perturbSS[name][i - 1] = experiment.perturbSS[name][0]

        # Write to file
        fn = join(self.resultsDir, self.name + '_equil_steadyStates.csv')
        _writeDataToCSV(fn, self.equilSS)

        #fn = join(self.resultsDir, self.name + '_perturb_steadyStates.csv')
        #_writeDataToCSV(fn, self.perturbSS)
        return

    def calcCrossConditionCorrelations(self):
        # Pair up molecules for correlations
        names = self.experiments[0].equilSS.dtype.names
        ceRNAs = [name for name in names if 'ceRNA' in name]
        miRNAs = [name for name in names if 'miRNA' in name]
        cePairs = list(combinations(ceRNAs, 2))
        miPairs = list(combinations(miRNAs, 2))

        # Create tables for results
        self.ceEquilCCs = np.zeros(len(cePairs), dtype=[('mol pair', 'a20'),
                                                ('r', 'f8'), ('p', 'f8')])
        #self.cePerturbCCs = np.zeros(len(cePairs), dtype=[('mol pair', 'a20'),
                                                #('r', 'f8'), ('p', 'f8')])
        self.miEquilCCs = np.zeros(len(miPairs), dtype=[('mol pair', 'a20'),
                                                ('r', 'f8'), ('p', 'f8')])
        #self.miPerturbCCs = np.zeros(len(miPairs), dtype=[('mol pair', 'a20'),
                                                #('r', 'f8'), ('p', 'f8')])

        # Do the calculations

        # ceRNAs
        # Equil steady states
        count = 0
        for pair in cePairs:
            mol1 = pair[0]
            mol2 = pair[1]
            mol1data = self.equilSS[mol1]
            mol2data = self.equilSS[mol2]

            (r, p) = pearsonr(mol1data, mol2data)

            if math.isnan(r):
                r = 0.0

            pairString = '(%s, %s)' % (mol1, mol2)
            self.ceEquilCCs['mol pair'][count] = pairString
            self.ceEquilCCs['r'][count] = r
            self.ceEquilCCs['p'][count] = p

            count += 1

        # Perturb steady states
        #count = 0
        #for pair in cePairs:
            #mol1 = pair[0]
            #mol2 = pair[1]
            #mol1data = self.perturbSS[mol1]
            #mol2data = self.perturbSS[mol2]

            #(r, p) = pearsonr(mol1data, mol2data)

            #if math.isnan(r):
                #r = 0.0

            #pairString = '(%s, %s)' % (mol1, mol2)
            #self.cePerturbCCs['mol pair'][count] = pairString
            #self.cePerturbCCs['r'][count] = r
            #self.cePerturbCCs['p'][count] = p

            #count += 1

        # miRNAs
        # Equil steady states
        count = 0
        for pair in miPairs:
            mol1 = pair[0]
            mol2 = pair[1]
            mol1data = self.equilSS[mol1]
            mol2data = self.equilSS[mol2]

            (r, p) = pearsonr(mol1data, mol2data)

            if math.isnan(r):
                r = 0.0

            pairString = '(%s, %s)' % (mol1, mol2)
            self.miEquilCCs['mol pair'][count] = pairString
            self.miEquilCCs['r'][count] = r
            self.miEquilCCs['p'][count] = p

            count += 1

        # Perturb steady states
        #count = 0
        #for pair in miPairs:
            #mol1 = pair[0]
            #mol2 = pair[1]
            #mol1data = self.perturbSS[mol1]
            #mol2data = self.perturbSS[mol2]

            #(r, p) = pearsonr(mol1data, mol2data)

            #if math.isnan(r):
                #r = 0.0

            #pairString = '(%s, %s)' % (mol1, mol2)
            #self.miPerturbCCs['mol pair'][count] = pairString
            #self.miPerturbCCs['r'][count] = r
            #self.miPerturbCCs['p'][count] = p

            #count += 1

        # Write results to files
        fn = join(self.resultsDir, self.name + '_ceRNA_equil_CCs.csv')
        _writeDataToCSV(fn, self.ceEquilCCs, frmt=('%20s', '%.18e', '%.18e'))

        #fn = join(self.resultsDir, self.name + '_ceRNA_perturb_CCs.csv')
        #_writeDataToCSV(fn, self.cePerturbCCs, frmt=('%20s', '%.18e', '%.18e'))

        fn = join(self.resultsDir, self.name + '_miRNA_equil_CCs.csv')
        _writeDataToCSV(fn, self.miEquilCCs, frmt=('%20s', '%.18e', '%.18e'))

        #fn = join(self.resultsDir, self.name + '_miRNA_perturb_CCs.csv')
        #_writeDataToCSV(fn, self.miPerturbCCs, frmt=('%20s', '%.18e', '%.18e'))

        # Make histograms
        fp = join(self.resultsDir, self.name + '_ceRNA_equil_CCs_hist')
        _histogram(fp, self.ceEquilCCs['r'], 'r')

        #fp = join(self.resultsDir, self.name + '_ceRNA_perturb_CCs_hist')
        #_histogram(fp, self.cePerturbCCs['r'], 'r')

        fp = join(self.resultsDir, self.name + '_miRNA_equil_CCs_hist')
        _histogram(fp, self.miEquilCCs['r'], 'r')

        #fp = join(self.resultsDir, self.name + '_miRNA_perturb_CCs_hist')
        #_histogram(fp, self.miPerturbCCs['r'], 'r')
        return

    def runAnalyses(self):
        self.collectSteadyStates()
        self.calcCrossConditionCorrelations()
        return

    def runAll(self):
        tStart = time.time()
        self.runExperiments()
        self.runAnalyses()
        tEnd = time.time()
        tElapsed = tEnd - tStart
        if _logging:
            self.logger.info('Time elapsed %.3f' % tElapsed)
        #print('Time elapsed %.3f' % tElapsed)
        return


class Population:
    def __init__(self, p, m, n, k, s, method, tEnd, outFreq, paramDict,
                 timestamp=None):

        if timestamp is None:
            self.timestamp = genTimeString()
        else:
            self.timestamp = timestamp

        self.m = m
        self.n = n
        self.k = k
        self.s = s
        self.p = p
        self.name = 'Pop@%d_%dx%dc%dx%d' % (p, m, n, k, s)
        self.createDirectories()

        if _logging:
            self.createLoggers()
            msg = 'Initializing population of size %d, with %d %dx%dx%d ensembles ' % (p, s, m, n, k)
            self.logger.info(msg)

        self.createEnsembles(m, n, k, s, method, tEnd, outFreq, paramDict)
        self.runAll()
        self.runAnalysis()
        return

    def createDirectories(self):
        self.rootDir = join(os.path.expanduser('~'),
                            'research/results/ceRNA/endysi/' + self.name)

        self.curRun = join(self.rootDir, self.timestamp)
        self.dataDir = join(self.curRun, 'data')
        self.resultsDir = join(self.curRun, 'results')
        makeDirs(self.dataDir)
        makeDirs(self.resultsDir)
        return

    def createEnsembles(self, m, n, k, s, method, tEnd, outFreq, paramDict):
        if _logging:
            self.logger.info('Creating ensembles')
        self.ensembles = []
        for i in range(self.p):
            self.ensembles.append(Ensemble(m, n, k, s, method, tEnd,
                                  outFreq, paramDict, baseDir=self.dataDir,
                                  timestamp='e%d' % (i + 1)))
        return

    def createLoggers(self):
        # create loggers and set levels
        self.debugger = logging.getLogger('debugger')
        self.debugger.setLevel(_debugLevel)
        self.logger = logging.getLogger('logger')
        self.logger.setLevel(logging.INFO)

        # create file handler which writes important events to a log file
        filename = join(self.curRun, 'pop_log.log')
        fh = logging.FileHandler(filename, mode='w')
        fh.setLevel(logging.INFO)

        # create console handler for writing debugging messages to the console
        ch = logging.StreamHandler()
        ch.setLevel(_debugLevel)

        # create formatter and add it to the handlers
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        # add the handlers to the logger
        self.logger.addHandler(fh)
        self.debugger.addHandler(ch)
        return

    def runAll(self):
        i = 1
        for e in self.ensembles:
            print('running ensemble %d' % i)
            if _logging:
                self.logger.info('Running simulations on ' + e.name)
            e.runAll()
            i += 1

        self.runAnalysis()
        return

    def runAnalysis(self):
        rVals = []
        for e in self.ensembles:
            fn = e.name + '_ceRNA_equil_CCs.csv'
            da = np.genfromtxt(join(e.resultsDir, fn), delimiter=';', names=True)
            rVals.extend(da['r'])

        fp = join(self.resultsDir, self.name + '_allCorrs')
        _histogram(fp, rVals, 'r')

        print('# of r vals = %d' % len(rVals))
        return


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', type=int, default=1,
                        help='The size of the population of ensembles')
    parser.add_argument('-m', type=int, help='The number of miRNAs')
    parser.add_argument('-n', type=int, help='The number of ceRNAs')
    parser.add_argument('-k', type=int, help='The degree of each ceRNA')
    parser.add_argument('-s', type=int,
                        help='The number of models in each ensemble')
    parser.add_argument('--method', type=str,
                        help='Simulation method: ode or ssa')
    args = parser.parse_args()

    paramDict = {}
    paramDict['vol'] = 2.0e-12
    paramDict['pT'] = (2.4e-3, 2.4e-1)
    paramDict['pR'] = (2.4e-3, 2.4e-1)
    paramDict['dT'] = (2.5e-05, 2.5e-03)
    paramDict['dR'] = (1e-05, 1e-03)
    paramDict['b'] = (1e-04, 1e-02)
    paramDict['u'] = (0, 1)
    paramDict['c'] = (7e-3, 7e-2)
    paramDict['a'] = (0.5, 0.5)

    maxHalfLife = 70000
    halfLifeMults = 2
    outFreq = 10
    nSamples = 2

    if args.method == 'ssa':
        outFreq = 100
        nSamples = 300

    tEnd = maxHalfLife * halfLifeMults * nSamples

    if args.p == 1:
        eds = Ensemble(args.m, args.n, args.k, args.s, args.method, tEnd,
                     outFreq, paramDict)
        eds.runAll()
    else:
        p = Population(args.p, args.m, args.n, args.k, args.s, args.method,
                       tEnd, outFreq, paramDict)
        p.runAll()
