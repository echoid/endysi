#!/usr/bin/env python2

# This script creates a 10x10 ceRNET with 10 links, and repeats experiments
# using different miRNA transcription rates, converting the NitzanParameters
# into the NitzanParametersReduced


import numpy as np
import endysi
import parameters as ps

m = n = k = 10
s = p = 100
outFreq = 2000
method = 'ode'

nSteps = 10

mins = np.logspace(np.log10(0.00024), np.log10(0.0024), num=nSteps)
maxs = np.logspace(np.log10(0.024), np.log10(0.24), num=nSteps)
minmax = zip(mins, maxs)

for i in range(nSteps):
    mm = minmax[i]
    params = ps.NitzanParametersCustom(mm[0], mm[1])
    endysi.makeRocketGoNow(m, n, k, s, p, outFreq, method, params)
