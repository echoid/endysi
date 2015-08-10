#!/usr/bin/env python2

# Module for interactive experimentation with ceRNETs

from __future__ import print_function, division
import os
import math
import time
import random
import numpy as np
import bngl
import endysi
from utilities import *


def initialize():
    createDirectories()


def createDirectories():
    rootDir = join(os.path.expanduser('~'),
                   'research/results/ceRNA/endysi/interactive')

    print(rootDir)


def newModel(m, n, k, method='ode'):
    pass


def setParameter(param, val):
    pass


def saveParameters(label='params'):
    pass


def setConcentration(mol, conc):
    pass


def saveConcentrations(label='concs'):
    pass


def writeConcentrations(suffix):
    pass


def writeParameters(suffix):
    pass


def run():
    running = True
    while running:
        # do some shit
        pass
