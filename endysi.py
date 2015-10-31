#!/usr/bin/env python
###############################################################################
# Endysi: ensemble dynamics simulator for ceRNA networks
###############################################################################

from __future__ import print_function, division
import os
import math
import time
import random
import socket
import logging
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import bngl
import parameters as params
from utilities import *

# Global debugging level
_debugLevel = logging.INFO  # use logging.INFO to turn off debugging

# Global logging switch
_logging = False

# Global dpi for plots
_dpi = 120

# Global setting for random seeding
_seeding = True

# global setting for alpha ranging
_rangingAlpha = False

_plottingTrajectories = False

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

colours10 = [purple, magenta, pink, red, yellow, green, blue,
           navy, aqua, orange]


# Global functions (mostly for analysis and plotting)
def _plotAllTrajectories(filename, data, mols, colours=None, ddpi=120):
    """Plot the trajectories of a set of molecules.  Assumes that the
    names in the mols argument are of the form 'mol_free', as in the
    output files from the BNG simulation.  Also assumes that the
    filename argument is the full absolute path for the file."""

    if len(mols) == 0:
        return

    if colours is None:
        colours = colours10

    count = 0
    fig = plt.figure()
    for mol in mols:
        plt.plot(data['time'], data[mol], colours[count],
                 label=mol.split('_')[0])
        count += 1

    plt.xlabel('Time (s)')
    plt.ylabel('Free molecules')
    plt.legend()
    plt.savefig(filename, dpi=_dpi, bbox_inches='tight')
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
    plt.savefig(filepath + '.png', dpi=_dpi, bbox_inches='tight')
    plt.close(fig)

    # Save the data
    histData = np.zeros(nBins, dtype=[('bin', 'f8'), ('count', 'f8')])
    np.copyto(histData['bin'], bins[:len(counts)])
    np.copyto(histData['count'], counts)
    _writeDataToCSV(filepath + '_histData.csv', histData)


def _scatterPlot(filename, x, y, coloured, xLabel, yLabel, size=30, title=None):
    if coloured:
        colour = colours10
    else:
        colour = 'b'

    fig = plt.figure()
    plt.scatter(x, y, s=size, c=colour, linewidth=0.1)
    if title is not None:
        fig.suptitle(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.savefig(filename, dpi=_dpi, bbox_inches='tight')
    plt.close(fig)


def _writeDataToCSV(filename, data, delim=';', comment='', frmt='%.18e'):
    head = ''
    for name in data.dtype.names:
        head += '%s;' % name
    head = head[:-1]

    np.savetxt(filename, data, delimiter=delim, header=head,
               comments=comment, fmt=frmt)


def _writeListToCSV(filename, data, name):
    with open(filename, 'w') as f:
        f.write('%s;\n' % name)
        for item in data:
            f.write('%f;\n' % item)


def _autocorrelate(X, maxDelta, outArray, freq=1, offset=0):
    """Calculate correlations of an array of values with itself.

    X is the input data array.
    maxDelta is the maximum delta (i.e., t - delta).
    outArray is the array in which to store resulting r values.
    freq is an optional output frequency factor.
    """

    for d in range(1, maxDelta + 1):
        delta = d * freq
        (r, p) = pearsonr(X[offset:-delta], X[delta + offset:])
        outArray[d - 1] = r

    return


### Class definitions ###
class Experiment:
    def __init__(self, parent, method, model, tEnd, outFreq, nSamples,
                 pTarget='ceRNA', pProp=5.0):

        self.model = model
        #self.simulator = bngl.BnglSimulator(model)
        self.perturbTarget = pTarget
        self.perturbProp = pProp
        self.method = method
        self.outFreq = outFreq
        self.nSamples = nSamples
        self.parent = parent

        #self.action = 'simulate({method=>"%(meth)s",suffix=>"%(suf)s",' + \
                      #'continue=>%(cnt)d,steady_state=>%(ss)d,' + \
                      #'t_end=>%(tend)d,n_steps=>%(nsteps)d,' + \
                      #'print_CDAT=>%(pcdat)d})'

        if method == 'ode':
            self.action = 'run_network -o ./%(name)s_%(suf)s -p cvode ' + \
                           '-a 1e-08 -r 1e-08 -c --cdat %(pcdat)d --fdat 0 ' + \
                           '-g ./%(name)s.net ./%(name)s.net %(outFreq)d ' + \
                           '%(tend)d'
        elif method == 'ssa':
            if _seeding:
                self.action = 'run_network -o ./%(name)s_%(suf)s -p ssa -h ' + \
                               '%(seed)d --cdat %(pcdat)d --fdat 0 -g ' + \
                               './%(name)s.net ./%(name)s.net ' + \
                               '%(outFreq)d %(tend)d'
            else:
                self.action = 'run_network -o ./%(name)s_%(suf)s -p ssa ' + \
                               '--cdat %(pcdat)d --fdat 0 -g ' + \
                               './%(name)s.net ./%(name)s.net ' + \
                               '%(outFreq)d %(tend)d'

        #nsteps = int(tEnd / outFreq)
        self.opts = {'suf': 'equil', 'tend': tEnd, 'outFreq': outFreq, 'cnt': 0,
                     'meth': method, 'pcdat': 0, 'ss': 0, 'name': model.name,
                     'seed': self.model.seed}

        if method == 'ode':
            self.opts['ss'] = 1
            self.opts['tend'] = 100000000000

        # fix tEnd business
        self.opts['tend'] = int(self.opts['tend'] / outFreq)

        return

    def purge(self):
        #del self.simulator
        del self.action
        del self.perturbProp
        del self.perturbTarget
        del self.opts
        self.model.purge()

    def run(self):
        os.chdir(self.model.home)

        # Equilibration phase
        self.opts['suf'] = 'equil'
        a = self.action % self.opts

        fn = self.model.filePath + '_sim_log.log'
        os.system(a + ' > ' + fn)

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
        fn = self.model.filePath + '_steadyStates.csv'
        _writeDataToCSV(fn, self.equilSS)

        return

    def calcWithinConditionCorrelations(self):
        # Pair up molecules for correlations
        cePairs = list(combinations(self.ceRNAs, 2))
        if self.parent.m > 1:
            miPairs = list(combinations(self.miRNAs, 2))

        # Create tables for results
        self.ceEquilWCCs = np.zeros(len(cePairs), dtype=[('mol pair', 'a20'),
                                    ('r', 'f8'), ('p', 'f8')])

        if self.parent.m > 1:
            self.miEquilWCCs = np.zeros(len(miPairs), dtype=[
                                        ('mol pair', 'a20'), ('r', 'f8'),
                                        ('p', 'f8')])

        sampFreq = int(140000 / self.outFreq)
        samplePoints = list(range(sampFreq, self.nSamples * sampFreq, sampFreq))

        # Create a directory for NaN correlation data
        self.nanDir = join(self.parent.resultsDir, 'nanData')
        makeDirs(self.nanDir)

        # ceRNAs
        count = 0
        nans = 0
        for pair in cePairs:
            mol1 = pair[0]
            mol2 = pair[1]
            mol1data = []
            mol2data = []

            # grab data at the sample frequency
            for i in samplePoints:
                mol1data.append(self.equilData[mol1][i])
                mol2data.append(self.equilData[mol2][i])

            (r, p) = pearsonr(mol1data, mol2data)

            if math.isnan(r):
                nans += 1
                r = 0.0
                #nanData = np.zeros(self.nSamples - 1, dtype=[(mol1, 'f8'),
                                                         #(mol2, 'f8')])
                #nanData[mol1] = mol1data
                #nanData[mol2] = mol2data
                #fn = join(self.nanDir, '%s_%s_%s.csv' % (self.model.name, mol1,
                                                         #mol2))
                #_writeDataToCSV(fn, nanData, frmt=('%.18e', '%.18e'))

            pairString = '(%s, %s)' % (mol1, mol2)
            self.ceEquilWCCs['mol pair'][count] = pairString
            self.ceEquilWCCs['r'][count] = r
            self.ceEquilWCCs['p'][count] = p

            count += 1

        # miRNAs
        count = 0
        nans = 0
        if self.parent.m > 1:
            for pair in miPairs:
                mol1 = pair[0]
                mol2 = pair[1]
                mol1data = []
                mol2data = []

                # grab data at the sample frequency
                for i in samplePoints:
                    mol1data.append(self.equilData[mol1][i])
                    mol2data.append(self.equilData[mol2][i])

                (r, p) = pearsonr(mol1data, mol2data)

                if math.isnan(r):
                    nans += 1
                    r = 0.0
                    #nanData = np.zeros(self.nSamples - 1, dtype=[(mol1, 'f8'),
                                            #(mol2, 'f8')])
                    #nanData[mol1] = mol1data
                    #nanData[mol2] = mol2data
                    #fn = join(self.nanDir, '%s_%s_%s.csv' % (self.model.name,
                                                             #mol1, mol2))
                    #_writeDataToCSV(fn, nanData, frmt=('%.18e', '%.18e'))

                pairString = '(%s, %s)' % (mol1, mol2)
                self.miEquilWCCs['mol pair'][count] = pairString
                self.miEquilWCCs['r'][count] = r
                self.miEquilWCCs['p'][count] = p

                count += 1

        # Write results to files
        fn = self.model.filePath + '_ceRNA_WCCs.csv'
        _writeDataToCSV(fn, self.ceEquilWCCs, frmt=('%20s', '%.18e', '%.18e'))

        if self.parent.m > 1:
            fn = self.model.filePath + '_miRNA_WCCs.csv'
            _writeDataToCSV(fn, self.miEquilWCCs, frmt=('%20s', '%.18e',
                                                        '%.18e'))

        return

    def calcAutocorrelations(self):
        maxHalfLife = 70000
        halfLifeMults = 2
        simTime = maxHalfLife * halfLifeMults * self.nSamples
        sampFreq = int((simTime / self.outFreq) / self.nSamples)
        maxDelta = 50
        startPoint = sampFreq

        # create array for r values
        dt = [(mol.split('_')[0], 'f8') for mol in self.ceRNAs]
        dt.insert(0, ('delta', 'i4'))
        rVals = np.zeros(maxDelta, dtype=dt)
        rVals['delta'] = np.arange(1, maxDelta + 1)

        # calculate autocorrs
        for mol in self.ceRNAs:
            molName = mol.split('_')[0]
            _autocorrelate(self.equilData[mol], maxDelta, rVals[molName],
                           freq=sampFreq, offset=startPoint)

        # write ACs to file
        head = 'delta'
        for mol in self.ceRNAs:
            head += ';%s' % mol.split('_')[0]
        fn = self.model.filePath + '_autocorrelations.csv'
        np.savetxt(fn, rVals, delimiter=';', header=head, comments='')
        return

    def plotTrajectories(self, ddpi=120):
        # Plot miRNAs
        fn = join(self.model.plotDir, self.model.name + '_miRNA_traj')
        _plotAllTrajectories(fn, self.equilData, self.miRNAs)

        #fn = join(self.model.plotDir, self.model.name + '_miRNA_perturb_traj')
        #_plotAllTrajectories(fn, self.perturbData, self.miRNAs)

        # Plot ceRNAs
        fn = join(self.model.plotDir, self.model.name + '_ceRNA_traj')
        _plotAllTrajectories(fn, self.equilData, self.ceRNAs)
        return

    def writeBindingPartners(self):
        fn = join(self.model.home, 'bindingPartners.txt')
        with open(fn, 'w') as f:
            f.write('Binding Partners:\n')
            f.write('ceRNAs:\n')
            for ceRNA in self.model.ceRNAs:
                p = [mol.name for mol in ceRNA.partners]
                f.write(ceRNA.name + ': ' + str(p) + '\n')

            f.write('\nmiRNAs:\n')
            for miRNA in self.model.miRNAs:
                p = [mol.name for mol in miRNA.partners]
                f.write(miRNA.name + ': ' + str(p) + '\n')

        return

    def gatherParams(self):
        pass

    def makeStochScatterPlots(self):
        if self.parent.n > 2:
            return

        sampFreq = int(140000 / self.outFreq)
        samplePoints = list(range(sampFreq, (self.nSamples + 1) * sampFreq,
                            sampFreq))

        dt = [('ceRNA1', 'f8'), ('ceRNA2', 'f8')]
        molData = np.zeros(len(samplePoints), dtype=dt)

        # grab data at the sample frequency
        count = 0
        for i in samplePoints:
            molData['ceRNA1'][count] = self.equilData['ceRNA1'][i]
            molData['ceRNA2'][count] = self.equilData['ceRNA2'][i]
            count += 1

        # Write results to files
        fn = self.model.filePath + '_ceRNA_samples.csv'
        _writeDataToCSV(fn, molData, frmt=('%.18e', '%.18e'))

        #t = 'ceRNA1 vs. ceRNA2'
        #fn = join(self.model.plotDir, self.model.name + '_ceRNA_scatter.png')
        #_scatterPlot(fn, molData['ceRNA1'], molData['ceRNA2'], False,
                     #'ceRNA1', 'ceRNA2', title=t)

    def makeScatterPlots(self):
        # gather params
        dt = [('steady state', 'f8'), ('trans', 'f8'), ('decay', 'f8'),
              ('halflife', 'f8')]
        ceParams = np.zeros(self.model.n, dtype=dt)
        miParams = np.zeros(self.model.m, dtype=dt)

        for ceRNA in self.model.ceRNAs:
            ceParams['trans'][ceRNA.num - 1] = ceRNA.prodRate
            ceParams['decay'][ceRNA.num - 1] = ceRNA.decayRate
            ceParams['halflife'][ceRNA.num - 1] = math.log(2) / ceRNA.decayRate
            ceParams['steady state'][ceRNA.num - 1] = \
                    self.equilSS[0][ceRNA.id - 1]

        for miRNA in self.model.miRNAs:
            miParams['trans'][miRNA.num - 1] = miRNA.prodRate
            miParams['decay'][miRNA.num - 1] = miRNA.decayRate
            miParams['halflife'][miRNA.num - 1] = math.log(2) / miRNA.decayRate
            miParams['steady state'][miRNA.num - 1] = \
                    self.equilSS[0][miRNA.id - 1]

        # save the data
        fn = self.model.filePath + '_ceRNA_ss_and_params.csv'
        _writeDataToCSV(fn, ceParams)

        fn = self.model.filePath + '_miRNA_ss_and_params.csv'
        _writeDataToCSV(fn, miParams)

        # make the plots
        # trans vs ceSS
        fn = join(self.model.plotDir, self.model.name + '_trans_vs_ceSS.png')
        t = 'ceRNAs: transcription rates vs steady states'
        _scatterPlot(fn, ceParams['trans'], ceParams['steady state'], True,
                     'transcription rate', 'ceRNA', title=t)

        # decay vs ceSS
        fn = join(self.model.plotDir, self.model.name + '_decay_vs_ceSS.png')
        t = 'ceRNAs: decay rates vs steady states'
        _scatterPlot(fn, ceParams['trans'], ceParams['steady state'], True,
                     'decay rate', 'ceRNA', title=t)

        # half-life vs ceSS
        fn = join(self.model.plotDir, self.model.name + '_halflife_vs_ceSS.png')
        t = 'ceRNAs: half-life vs steady states'
        _scatterPlot(fn, ceParams['halflife'], ceParams['steady state'], True,
                     'half-life', 'ceRNA', title=t)

        # miRNAs
        # trans vs miSS
        fn = join(self.model.plotDir, self.model.name + '_trans_vs_miSS.png')
        t = 'miRNAs: transcription rates vs steady states'
        _scatterPlot(fn, miParams['trans'], miParams['steady state'], True,
                     'transcription rate', 'miRNA', title=t)

        # decay vs miSS
        fn = join(self.model.plotDir, self.model.name + '_decay_vs_miSS.png')
        t = 'miRNAs: decay rates vs steady states'
        _scatterPlot(fn, miParams['trans'], miParams['steady state'], True,
                     'decay rate', 'miRNA', title=t)

        # half-life vs miSS
        fn = join(self.model.plotDir, self.model.name + '_halflife_vs_miSS.png')
        t = 'miRNAs: half-life vs steady states'
        _scatterPlot(fn, miParams['halflife'], miParams['steady state'], True,
                     'half-life', 'miRNA', title=t)

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
        fn = join(self.model.plotDir, self.model.name + '_foldPlots.png')
        plt.savefig(fn, dpi=ddpi, bbox_inches='tight')
        plt.close(fig)
        return

    def runAnalyses(self):
        self.loadData()
        self.calcSteadyStates()
        if self.method == 'ssa':
            self.calcWithinConditionCorrelations()
            self.calcAutocorrelations()
            self.makeStochScatterPlots()

        if _plottingTrajectories:
            self.plotTrajectories()

        #self.makeScatterPlots()
        self.writeBindingPartners()
        #self.makeFoldPlots()
        #self.deleteData()
        return

    def deleteData(self):
        del self.equilData
        return


class Ensemble:
    def __init__(self, m, n, k, size, method, tEnd, outFreq, nSamples,
                 randParams, timestamp=None, baseDir=None, seedScale=None,
                 alpha=None, linearSampling=0, name=None):

        if timestamp is None:
            self.timestamp = genTimeString()
        else:
            self.timestamp = timestamp

        if name is None:
            self.name = 'ceRNET_%dx%dc%dx%d' % (m, n, k, size)
        else:
            self.name = name

        self.m = m
        self.n = n
        self.k = k
        self.size = size
        self.method = method
        self.tEnd = tEnd
        self.outFreq = outFreq
        self.models = []
        self.experiments = []
        self.nSamples = nSamples
        self.linearSampling = bool(linearSampling)
        self.seedScale = seedScale
        self.randParams = randParams
        self.alpha = alpha

        self.createDirectories(baseDir)

        if _logging:
            self.createLoggers()
            msg = 'Initializing system to run {0} networks with M={1}, ' + \
                   'N={2}, K={3}'
            self.logger.info(msg.format(size, m, n, k))
            self.logger.info('Parameter ranges are as follows:')
            self.logger.info(str(randParams))

            self.logger.info('Creating models...')

        self.writeRunInfo()
        return

    def writeRunInfo(self):
        with open(join(self.resultsDir, 'runInfo'), 'w') as riFile:
            riFile.write('Ensemble of %d parameter sets\n' % self.size)
            riFile.write('m = %d\n' % self.m)
            riFile.write('n = %d\n' % self.n)
            riFile.write('k = %d\n' % self.k)
            riFile.write('method: %s\n' % self.method)
            riFile.write('Parameters: %s\n' % self.randParams.name)

        with open(join(self.resultsDir, 'alpha'), 'w') as aFile:
            aFile.write('alpha=' + str(self.alpha))

        return

    def purge(self):
        #del self.randParams
        del self.tEnd
        del self.linearSampling
        del self.m
        del self.n
        del self.k
        del self.size
        del self.method
        del self.nSamples

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
                                'research/sims/ceRNA/endysi/')
        else:
            self.rootDir = baseDir

        #self.runDir = join(self.rootDir, self.name)
        bResultsDir = join(self.rootDir, 'results')
        bDataDir = join(self.rootDir, 'data')
        self.rRunDir = join(bResultsDir, self.name)
        self.dRunDir = join(bDataDir, self.name)
        self.resultsDir = join(self.rRunDir, self.timestamp)
        self.dataDir = join(self.dRunDir, self.timestamp)
        makeDirs(self.dataDir)
        makeDirs(self.resultsDir)
        return

    def createModels(self, randParams):
        for i in range(1, self.size + 1):
            dDir = join(self.dataDir, 'model%d' % i)
            model = bngl.CernetModel(dDir, self.m, self.n, self.k, i,
                                     randParams, seed=None)

            self.models.append(model)
        return

    def createExperiments(self, method, tEnd, outFreq):
        for model in self.models:
            e = Experiment(method, model, tEnd, outFreq)
            self.experiments.append(e)

    def runExperiments(self):
        for experiment in self.experiments:
            if _logging:
                m = 'Running simulations on %s' % experiment.model.name
                self.logger.info(m)

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

        for experiment in self.experiments:
            i = experiment.model.index
            self.equilSS['model'][i - 1] = i
            #self.perturbSS['model'][i - 1] = i
            for name in names:
                self.equilSS[name][i - 1] = experiment.equilSS[name][0]

        # Write to file
        fn = join(self.resultsDir, self.name + '_steadyStates.csv')
        _writeDataToCSV(fn, self.equilSS)

        return

    def calcBalance(self):
        # Calc RNA means
        names = self.experiments[0].equilSS.dtype.names
        dt = [('model', 'i4'), ('ceMean', 'f8'), ('miMean', 'f8'),
              ('mean_ceCorr', 'f8'), ('mean_miCorr', 'f8'), ('ratio', 'f8'),
              ('ceRxceMean', 'f8'), ('miRxmiMean', 'f8'), ('integ', 'f8')]
        self.nSSandR = np.zeros(self.size, dtype=dt)

        for experiment in self.experiments:
            # calc norm'd avg steady states
            i = experiment.model.index
            self.nSSandR['model'][i - 1] = i
            miTot = 0
            ceTot = 0
            ceCount = 0.0
            miCount = 0.0

            for name in names:
                if 'mi' in name:
                    miTot += experiment.equilSS[name][0]
                    miCount += 1.0
                elif 'ce' in name:
                    ceTot += experiment.equilSS[name][0]
                    ceCount += 1.0

            self.nSSandR['miMean'][i - 1] = miTot / miCount
            self.nSSandR['ceMean'][i - 1] = ceTot / ceCount

            # add correlations
            if self.method == 'ode':
                self.nSSandR['mean_ceCorr'][i - 1] = np.mean(
                                self.ceEquilCCs['r'])
                if self.m > 1:
                    self.nSSandR['mean_miCorr'][i - 1] = np.mean(
                                self.miEquilCCs['r'])

            elif self.method == 'ssa':
                self.nSSandR['mean_ceCorr'][i - 1] = np.mean(self.ceWCCCs)
                if self.m > 1:
                    self.nSSandR['mean_miCorr'][i - 1] = np.mean(self.miWCCCs)

            # calc other measures
            ceR = self.nSSandR['mean_ceCorr'][i - 1]
            miR = self.nSSandR['mean_miCorr'][i - 1]
            ceMean = self.nSSandR['ceMean'][i - 1]
            miMean = self.nSSandR['miMean'][i - 1]

            self.nSSandR['ratio'][i - 1] = ceMean / miMean
            self.nSSandR['ceRxceMean'][i - 1] = ca = ceR * ceMean
            self.nSSandR['miRxmiMean'][i - 1] = ma = miR * miMean
            self.nSSandR['integ'][i - 1] = ca / ma

        # write data to file
        fn = join(self.resultsDir, self.name + '_nSSandRs.csv')
        _writeDataToCSV(fn, self.nSSandR)

        # Only plot if this is not part of a population run
        if self.method == 'ssa':
            fig = plt.figure()
            plt.scatter(self.nSSandR['ceMean'], self.nSSandR['mean_ceCorr'],
                     color='b', label='ceRNA')

            if self.m > 1:
                plt.scatter(self.nSSandR['miMean'],
                        self.nSSandR['mean_miCorr'], color='r', label='miRNA')

            cType = ''
            if self.method == 'ode':
                cType = 'CCC'
            else:
                cType = 'WCC'

            t = 'Normalized Average RNA vs Average %s Correlation' % cType
            fig.suptitle(t)
            plt.xlabel('Average RNA steady state')

            yLabel = 'Average ' + cType
            plt.ylabel(yLabel)

            # Save it
            fn = join(self.resultsDir, self.name + '_nSSand%s.png' % cType)
            plt.savefig(fn, dpi=_dpi, bbox_inches='tight')
            plt.close(fig)

        # calc overall mean
        self.ceMean = np.mean(self.nSSandR['ceMean'])
        self.miMean = np.mean(self.nSSandR['miMean'])

        return

    def calcMolSums(self):
        names = self.experiments[0].equilSS.dtype.names
        dt = [('model', 'i4'), ('miTot', 'f8'), ('ceTot', 'f8')]
        equilSums = np.zeros(self.size, dtype=dt)

        for experiment in self.experiments:
            i = experiment.model.index
            equilSums['model'][i - 1] = i
            miTot = 0
            ceTot = 0
            for name in names:
                if 'mi' in name:
                    miTot += experiment.equilSS[name][0]
                elif 'ce' in name:
                    ceTot += experiment.equilSS[name][0]

            equilSums['miTot'][i - 1] = miTot
            equilSums['ceTot'][i - 1] = ceTot

        # normalize
        normSums = np.copy(equilSums)
        maxR = max(max(equilSums['miTot']), max(equilSums['ceTot']))

        for i in range(self.size):
            ovs = equilSums['miTot'][i]
            ovr = equilSums['ceTot'][i]
            normSums['miTot'][i] = ovs / maxR
            normSums['ceTot'][i] = ovr / maxR

        # Write to files
        fn = join(self.resultsDir, self.name + '_sums.csv')
        _writeDataToCSV(fn, equilSums)

        fn = join(self.resultsDir, self.name + '_normSums.csv')
        _writeDataToCSV(fn, normSums)

    def calcWithinConditionCorrelations(self):
        # i.e., stochastic fluctuation correlations
        ceRVals = []
        miRVals = None
        if self.m > 1:
            miRVals = []

        for e in self.experiments:
            #print(type(e.ceEquilWCCs['r'][0]))
            #if isinstance(e.ceEquilWCCs['r'][0], np.ndarray):
            ceRVals.extend(e.ceEquilWCCs['r'])
            #else:
                #ceRVals.append(e.ceEquilWCCs['r'])

            if self.m > 1:
                #if isinstance(e.miEquilWCCs['r'][0], np.ndarray):
                miRVals.extend(e.miEquilWCCs['r'])
                #else:
                    #miRVals.append(e.miEquilWCCs['r'])

        # plot and save
        # Only plot it this is not part of a population
        if self.method == 'ssa':
            fp = join(self.resultsDir, self.name + '_ceWCCs')
            _histogram(fp, ceRVals, 'r')

        fn = join(self.resultsDir, self.name + '_ceWCCs.csv')
        _writeListToCSV(fn, ceRVals, 'r')

        if self.m > 1:
            if self.method == 'ssa':
                fp = join(self.resultsDir, self.name + '_miWCCs')
                _histogram(fp, miRVals, 'r')

            fn = join(self.resultsDir, self.name + '_miWCCs.csv')
            _writeListToCSV(fn, miRVals, 'r')

        # store rVals temporarily
        self.ceWCCCs = np.array(ceRVals)
        if self.m > 1:
            self.miWCCCs = np.array(miRVals)

        return

    def calcAutocorrelations(self):
        maxDelta = 50
        allRvals = []
        for e in self.experiments:
            fn = e.model.filePath + '_autocorrelations.csv'
            da = np.genfromtxt(fn, delimiter=';', names=True)
            allRvals.append(da)

        avgRvals = np.zeros(maxDelta, dtype=allRvals[0].dtype)
        avgRvals['delta'] = np.arange(1, maxDelta + 1)

        mols = [name for name in avgRvals.dtype.names if name != 'delta']
        for mol in mols:
            molData = [rVals[mol] for rVals in allRvals]
            for i in range(maxDelta):
                rAvg = sum([rVals[i] for rVals in molData]) / len(molData)
                avgRvals[mol][i] = rAvg

        # write avgs to file
        head = 'delta'
        for mol in mols:
            head += ';%s' % mol.split('_')[0]
        fn = join(self.resultsDir, '%s_avgACs.csv' % self.name)
        np.savetxt(fn, avgRvals, delimiter=';', header=head, comments='')

        # plot ACs
        fig = plt.figure()
        for mol in mols:
            plt.plot(avgRvals['delta'], avgRvals[mol], label=mol)
        plt.xlabel('delta')
        plt.ylabel('r')
        plt.legend()
        fn = join(self.resultsDir, '%s_avgACs.png' % self.name)
        plt.savefig(fn, dpi=120, bbox_inches='tight')
        plt.close(fig)

        return

    def calcCrossConditionCorrelations(self):
        # Pair up molecules for correlations
        names = self.experiments[0].equilSS.dtype.names
        ceRNAs = [name for name in names if 'ceRNA' in name]
        miRNAs = [name for name in names if 'miRNA' in name]
        cePairs = list(combinations(ceRNAs, 2))
        if self.m >= 2:
            miPairs = list(combinations(miRNAs, 2))

        # Create tables for results
        dt = [('mol pair', 'a20'), ('r', 'f8'), ('p', 'f8')]
        self.ceEquilCCs = np.zeros(len(cePairs), dtype=dt)
        if self.m >= 2:
            self.miEquilCCs = np.zeros(len(miPairs), dtype=dt)

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

        # Write results to files
        fn = join(self.resultsDir, self.name + '_ceRNA_CCs.csv')
        _writeDataToCSV(fn, self.ceEquilCCs, frmt=('%20s', '%.18e', '%.18e'))

        # miRNAs
        # Equil steady states
        if self.m >= 2:
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

            # Write results to files
            fn = join(self.resultsDir, self.name + '_miRNA_CCs.csv')
            _writeDataToCSV(fn, self.miEquilCCs, frmt=('%20s', '%.18e',
                                                       '%.18e'))

        return

    def runAnalyses(self):
        self.collectSteadyStates()
        #self.calcMolSums()

        if self.method == 'ode':
            self.calcCrossConditionCorrelations()
        else:
            self.calcWithinConditionCorrelations()
            self.calcAutocorrelations()

        self.calcBalance()

        return

    def runAll(self):
        # create pRange for fixed param runs
        pRange = None
        if self.randParams.isFixed():
            #pRange = np.linspace(self.randParams.pMin, self.randParams.pMax,
                                 #num=self.size)
            pRange = np.logspace(np.log10(self.randParams.pMin),
                                 np.log10(self.randParams.pMax),
                                 num=self.size)

            # save the range to file
            s = self.name
            if self.randParams.isDirect():
                s += '_%s_range' % self.randParams.pName
            else:
                s += '_range'
            fn = join(self.resultsDir, s)
            np.savetxt(fn, pRange)
            self.pRange = pRange

        elif self.randParams.isSpecific():
            pRange = self.randParams.ps
            self.pRange = pRange

        tStart = time.time()
        for i in range(1, self.size + 1):
            dDir = join(self.dataDir, 'model%d' % i)

            s = None
            if _seeding:
                s = i
                if self.seedScale is not None:
                    s *= self.seedScale

            if pRange is not None:
                if self.randParams.isDirect():
                    self.randParams.pVal = pRange[i - 1]
                elif (self.randParams.isFixed()
                          and not self.randParams.isDirect()):
                    self.randParams.set(self.randParams.pRanged, pRange[i - 1])
                elif self.randParams.isSpecific():
                    self.randParams.set(self.randParams.pSet, pRange[i - 1])

            model = bngl.CernetModel(dDir, self.m, self.n, self.k, i,
                                     self.randParams, alpha=self.alpha,
                                     seed=s, linearSampling=self.linearSampling)

            e = Experiment(self, self.method, model, self.tEnd, self.outFreq,
                           self.nSamples)

            e.run()
            e.runAnalyses()
            e.deleteData()
            e.purge()
            self.experiments.append(e)

        self.runAnalyses()
        tEnd = time.time()
        tElapsed = tEnd - tStart
        if _logging:
            self.logger.info('Time elapsed %.3f' % tElapsed)

        return


class Population:
    def __init__(self, p, m, n, k, s, method, tEnd, outFreq, randParams,
                 fixedParams=None, timestamp=None, baseDir=None,
                 rangingAlpha=False, linearSampling=1):

        if timestamp is None:
            self.timestamp = genTimeString()
        else:
            self.timestamp = timestamp

        self.m = m
        self.n = n
        self.k = k
        self.s = s
        self.p = p
        self.method = method
        self.tEnd = tEnd
        self.outFreq = outFreq
        self.randParams = randParams
        self.fixedParams = fixedParams
        self.name = 'Pop@%d_%dx%dc%dx%d' % (p, m, n, k, s)
        self.ensembles = []
        self.createDirectories(baseDir)
        self.linearSampling = bool(linearSampling)
        self.rangingAlpha = rangingAlpha

        if _logging:
            self.createLoggers()
            msg = 'Initializing population of size %d, with %d ' + \
                   '%dx%dx%d ensembles ' % (p, s, m, n, k)
            self.logger.info(msg)
        self.writeRunInfo()

        return

    def writeRunInfo(self):
        with open(join(self.resultsDir, 'runInfo'), 'w') as riFile:
            riFile.write('Population of %d ensembles, each of size %d\n' %
                         (self.p, self.s))
            riFile.write('m = %d\n' % self.m)
            riFile.write('n = %d\n' % self.n)
            riFile.write('k = %d\n' % self.k)
            riFile.write('method: %s\n' % self.method)
            riFile.write('rangingAlpha: ' + str(self.rangingAlpha) + '\n')
            riFile.write('Parameter ranges: %s\n' % self.randParams.name)
        return

    def createDirectories(self, baseDir):
        if baseDir is None:
            self.rootDir = join(os.path.expanduser('~'),
                            'research/sims/ceRNA/endysi/')
        else:
            self.rootDir = join(baseDir, self.name)

        bResultsDir = join(self.rootDir, 'results')
        bDataDir = join(self.rootDir, 'data')
        rRunDir = join(bResultsDir, self.name)
        dRunDir = join(bDataDir, self.name)
        self.resultsDir = join(rRunDir, self.timestamp)
        self.dataDir = join(dRunDir, self.timestamp)
        makeDirs(self.dataDir)
        makeDirs(self.resultsDir)
        return

    def createEnsembles(self, m, n, k, s, method, tEnd, outFreq, randParams):
        if _logging:
            self.logger.info('Creating ensembles')
        self.ensembles = []
        for i in range(self.p):
            self.ensembles.append(Ensemble(m, n, k, s, method, tEnd,
                                  outFreq, randParams, baseDir=self.dataDir,
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
        tStart = time.time()
        offset = 0
        alphas = None
        if self.rangingAlpha:
            alphas = np.linspace(0.01, 1, self.p)

        for i in range(self.p):
            if _logging:
                self.logger.info('Running simulations on ' + e.name)

            sScale = None
            if _seeding:
                sScale = 999 + offset
                offset += 3

            a = None
            if self.rangingAlpha:
                a = alphas[i]

            e = Ensemble(self.m, self.n, self.k, self.s, self.method, self.tEnd,
                         self.outFreq, 1, self.randParams,
                         baseDir=self.dataDir, timestamp='e%d' % (i + 1),
                         seedScale=sScale, linearSampling=self.linearSampling,
                         alpha=a)
            e.runAll()
            e.purge()
            self.ensembles.append(e)

        self.runAnalysis()

        tEnd = time.time()
        tElapsed = tEnd - tStart
        print('Time elapsed: %f' % tElapsed)
        return

    def alphaAnalysis(self):
        dt = [('alpha', 'f8'), ('normR', 'f8'), ('normS', 'f8'), ('rAvg', 'f8')]
        da = np.zeros(self.p, dtype=dt)
        dao = np.zeros(self.p, dtype=dt)

        i = 0
        for e in self.ensembles:
            fn = e.name + '_ceRNA_CCs.csv'
            ccd = np.genfromtxt(join(e.resultsDir, fn), delimiter=';',
                               names=True)
            fs = e.name + '_steadyStates.csv'
            ssd = np.genfromtxt(join(e.resultsDir, fs), delimiter=';',
                                names=True)

            alpha = 0
            with open(join(e.curRun, 'alpha'), 'r') as af:
                alpha = float(af.readline().split('=')[1])

            da['alpha'][i] = alpha
            da['rAvg'][i] = np.mean(ccd['r'])
            dao['alpha'][i] = alpha
            dao['rAvg'][i] = np.mean(ccd['r'])

            totR = 0
            totS = 0
            for name in ssd.dtype.names:
                if 'miRNA' in name:
                    totS += sum(ssd[name])
                elif 'ceRNA' in name:
                    totR += sum(ssd[name])

            # Individual normalization
            denom = totR + totS
            normR = totR / denom
            normS = totS / denom
            da['normR'][i] = normR
            da['normS'][i] = normS

            dao['normR'][i] = totR
            dao['normS'][i] = totS
            i += 1

        ## Global normalization
        maxR = max(max(dao['normR']), max(dao['normS']))

        for i in range(self.p):
            ovs = dao['normS'][i]
            ovr = dao['normR'][i]
            dao['normS'][i] = ovs / maxR
            dao['normR'][i] = ovr / maxR

        # write data to file
        fn = join(self.resultsDir, self.name + '_normSums_individual.csv')
        _writeDataToCSV(fn, da)

        fn = join(self.resultsDir, self.name + '_normSums_global.csv')
        _writeDataToCSV(fn, dao)

        # plot normalized RNA levels and R against alpha
        fig = plt.figure()
        plt.plot(da['alpha'], da['normR'], label='normalized ceRNA')
        plt.plot(da['alpha'], da['normS'], label='normalized miRNA')
        plt.plot(da['alpha'], da['rAvg'], label='average r')
        plt.xlabel('alpha')
        plt.legend()
        fn = join(self.resultsDir, self.name + '_normSumsVsAlpha_indiv.png')
        plt.savefig(fn, dpi=_dpi, bbox_inches='tight')
        plt.close(fig)

        fig = plt.figure()
        plt.plot(dao['alpha'], dao['normR'], label='normalized ceRNA')
        plt.plot(dao['alpha'], dao['normS'], label='normalized miRNA')
        plt.plot(dao['alpha'], dao['rAvg'], label='average r')
        plt.xlabel('alpha')
        plt.legend()
        fn = join(self.resultsDir, self.name + '_normSumsVsAlpha_global.png')
        plt.savefig(fn, dpi=_dpi, bbox_inches='tight')
        plt.close(fig)

    def calcBalance(self):
        # Create a table for results
        dt = [('ens', 'i4'), ('norm avg ceRNA', 'f8'), ('norm avg miRNA', 'f8'),
              ('avg ceCCC', 'f8'), ('avg miCCC', 'f8'), ('ratio', 'f8'),
              ('ceRxceMean', 'f8'), ('miRxmiMean', 'f8'), ('integ', 'f8')]
        self.nSSandR = np.zeros(self.p, dtype=dt)

        # Gather means from the ensembles
        ceMeans = [e.ceMean for e in self.ensembles]
        miMeans = [e.miMean for e in self.ensembles]

        # Normalize and store
        for i in range(self.p):
            # Norm'd means
            self.nSSandR['ens'][i] = i + 1
            denom = ceMeans[i] + miMeans[i]
            self.nSSandR['norm avg ceRNA'][i] = ceMeans[i] / denom
            self.nSSandR['norm avg miRNA'][i] = miMeans[i] / denom

            # Mean correlation
            e = self.ensembles[i]
            self.nSSandR['avg ceCCC'][i] = np.mean(e.ceEquilCCs['r'])
            if self.m > 1:
                self.nSSandR['avg miCCC'][i] = np.mean(e.miEquilCCs['r'])

            ceMean = self.nSSandR['norm avg ceRNA'][i - 1]
            miMean = self.nSSandR['norm avg miRNA'][i - 1]
            ceR = self.nSSandR['avg ceCCC'][i - 1]
            miR = self.nSSandR['avg miCCC'][i - 1]

            self.nSSandR['ratio'][i - 1] = ceMean / miMean
            self.nSSandR['ceRxceMean'][i - 1] = ca = ceR * ceMean
            self.nSSandR['miRxmiMean'][i - 1] = ma = miR * miMean
            self.nSSandR['integ'][i - 1] = ca / ma

        # Save to file
        fn = join(self.resultsDir, self.name + '_normSSandR.csv')
        _writeDataToCSV(fn, self.nSSandR)

        # Plot it
        fig = plt.figure()
        plt.scatter(self.nSSandR['norm avg ceRNA'], self.nSSandR['avg ceCCC'],
                 color='b', label='ceRNA')
        plt.scatter(self.nSSandR['norm avg miRNA'], self.nSSandR['avg miCCC'],
                 color='r', label='miRNA')

        t = 'Normalized Average RNA vs Average cross-conditional Correlation'
        fig.suptitle(t)
        plt.xlabel('Average RNA steady state')
        plt.ylabel('Average CCC')

        # Save it
        fn = join(self.resultsDir, 'nSSandCCC.png')
        plt.savefig(fn, dpi=_dpi, bbox_inches='tight')
        plt.close(fig)

    def collectCrossConditionCorrelations(self):
        ceCCs = []
        miCCs = []

        for e in self.ensembles:
            fn = join(e.resultsDir, e.name + '_ceRNA_CCs.csv')
            ceda = np.genfromtxt(fn, delimiter=';', names=True)

            if self.n <= 2:
                ceCCs.append(ceda['r'])
            else:
                ceCCs.extend(ceda['r'])

        fp = join(self.resultsDir, self.name + '_ceCCCs')
        _histogram(fp, ceCCs, 'r')

        fn = join(self.resultsDir, self.name + '_ceCCCs.csv')
        _writeListToCSV(fn, ceCCs, 'r')

        if self.m >= 2:
            for e in self.ensembles:
                fn = join(e.resultsDir, e.name + '_miRNA_CCs.csv')
                mida = np.genfromtxt(fn, delimiter=';', names=True)

                if self.m >= 2:
                    if self.m == 2:
                        miCCs.append(mida['r'])
                    else:
                        miCCs.extend(mida['r'])

            fp = join(self.resultsDir, self.name + '_miCCCs')
            _histogram(fp, miCCs, 'r')

            fn = join(self.resultsDir, self.name + '_miCCCs.csv')
            _writeListToCSV(fn, miCCs, 'r')

        return

    def runAnalysis(self):
        self.collectCrossConditionCorrelations()
        self.calcBalance()

        if self.rangingAlpha:
            self.alphaAnalysis()

        return


def makeRocketGoNow(m, n, k, s, p, outFreq, method, randParams, linSamp=False,
                    rangeAlpha=False):

    maxHalfLife = 70000
    halfLifeMults = 2
    nSamples = 1

    if method == 'ssa':
        nSamples = 1000

    tEnd = maxHalfLife * halfLifeMults * nSamples

    # Check if we're running on the lab cluster
    baseDir = None
    host = socket.gethostname()
    if host == 'crick':
        baseDir = '/ohri/projects/perkins/mattm/ceRNA/endysi'
    elif host == 'ogic':
        baseDir = '/data/matt/ceRNA/endysi'

    if p == 1:
        eds = Ensemble(m, n, k, s, method, tEnd, outFreq, nSamples,
                       randParams, linearSampling=linSamp,
                       baseDir=baseDir)
        eds.runAll()
    else:
        pop = Population(p, m, n, k, s, method, tEnd, outFreq, randParams,
                       linearSampling=linSamp, baseDir=baseDir,
                       rangingAlpha=rangeAlpha)
        pop.runAll()

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
    parser.add_argument('--linear', type=int, default=0,
                        help='Whether to sample over linear or log space' +
                        'Default is log sampling, which works better')
    parser.add_argument('--alpha', type=int, default=0,
                        help='Range alpha parameter; default off (0)')
    parser.add_argument('-o', type=int, default=2000, help='Output frequency')
    parser.add_argument('-t', type=str, default='',
                        help='timestamp for restart')

    args = parser.parse_args()

    randParams = params.NitzanParametersReduced()
    makeRocketGoNow(args.m, args.n, args.k, args.s, args.p, args.o,
                    args.method, randParams, linSamp=args.linear,
                    rangeAlpha=args.alpha)
