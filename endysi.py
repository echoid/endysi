#!/usr/bin/env python

import os
import sys
import math
import time
import glob
import random
import logging
from itertools import combinations
import pexpect
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import bngl
from bngci import BNGCI
from utilities import *

# global debugging level:
_debugLevel = logging.DEBUG # use logging.INFO to turn off debugging


### Class definitions ###
class Simulation:
    def __init__(self, model, simulator, tEnd, outFreq):
        self.model = model
        self.sim = simulator
        


class Experiment:
    pass


class Endysi:
    def __init__(self, m, n, k, size, timestamp=None):
        if timestamp is None:
            self.timestamp = genTimeString()
        else:
            self.timestamp = timestamp
        
        self.name = 'ceRNET_%dx%d_%d' % (m, n, k)
        self.m = m
        self.n = n
        self.k = k
        self.size = size
        
        createDirectories()
        createLoggers()
        self.logger.info('Initializing Endysi to run {0} networks with M={1}, N={2}, K={3}'.format(size, m, n, k))
        
    
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
    
    def createDirectories(self):
        self.rootDir = join(os.path.expanduser('~'), 'research/results/ceRNA/' + self.name)
        self.curRun = join(self.rootDir, self.timestamp)
        self.dataDir = join(self.curRun, 'data')
        self.resultsDir = join(self.curRun, 'results')
        makeDirs(self.dataDir)
        makeDirs(self.resultsDir)
    
    def createModels(self):
        for i in range(1, self.size+1):
            





