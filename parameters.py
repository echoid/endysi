#!/usr/bin/env python2

###############################################################################
# Classes for parameter ranges
###############################################################################
from __future__ import print_function


class ParameterRange(object):
    def __init__(self, name, fixed=False, direct=False, specific=False):
        self.name = name
        self.params = {}
        self.fixed = fixed
        self.direct = direct
        self.specific = specific

    def get(self, param):
        if param in self.params:
            return self.params[param]
        else:
            print('Parameter %s not found!' % param)
            return None

    def getMin(self, param):
        p = self.get(param)
        if p is None:
            return None
        else:
            return p[0]

    def getMax(self, param):
        p = self.get(param)
        if p is None:
            return p

        if len(p) == 2:
            return p[1]
        else:
            print('Parameter %s is fixed; no max val available!' % param)
            return None

    def set(self, param, val):
        if isinstance(val, tuple):
            self.params[param] = val
        else:
            self.params[param] = (val,)

    def isFixed(self):
        return self.fixed

    def isDirect(self):
        return self.direct

    def isSpecific(self):
        return self.specific


class NitzanParameters(ParameterRange):
    """This class represents the parameter ranges for ceRNA networks from
    Nitzan et al., 2014.

    Values are presented in Table S1, column Figure 5.  The only change is that
    parameter u_c has been changed from 0 - 1 in the original to 1e-04 - 1e-02.
    """

    def __init__(self):
        super(NitzanParameters, self).__init__('NitzanParameters')
        self.params = {'pR': (2.4e-03, 2.4e-01), 'pS': (2.4e-03, 2.4e-01),
                       'dR': (1e-05, 1e-03), 'dS': (2.5e-05, 2.5e-03),
                       'b': (1e-04, 1e-02), 'u': (1e-04, 1e-02),
                       'c': (7e-03, 7e-02), 'a': (0.5, 0.5)}


class NitzanParametersReduced(NitzanParameters):
    """This class represents the parameter ranges from Nitzan et al., 2014.
    The range of miRNA transcription rates (pS) has been modified in an attempt
    to keep miRNA and ceRNA levels more balanced.

    Specifically, the original min and max values have been reduced by an order
    of magnitute.
    """

    def __init__(self):
        super(NitzanParametersReduced, self).__init__()
        self.set('pS', (2.4e-05, 2.4e-02))
        self.name = 'NitzanParametersReduced'


class NitzanParametersCustom(NitzanParametersReduced):

    def __init__(self, sMin, sMax):
        super(NitzanParametersCustom, self).__init__()
        self.set('pS', (sMin, sMax))
        self.name = 'NitzanParametersCustom_%f--%f' % (sMin, sMax)


class NitzanParametersExpanded(NitzanParameters):
    """This class represents an expanded parameter range based on Nitzan et al.,
    2014.

    The original ranges have been expanded by two orders of magnitude, extending
    it in both directions (min and max)
    """

    def __init__(self):
        super(NitzanParametersExpanded, self).__init__()
        self.params = {'pR': (2.4e-04, 2.4), 'pS': (2.4e-04, 2.4),
                       'dR': (1e-06, 1e-02), 'dS': (2.5e-06, 2.5e-02),
                       'b': (1e-05, 1e-01), 'u': (1e-05, 1e-01),
                       'c': (7e-04, 7e-01), 'a': (0.5, 0.5)}
        self.name = 'NitzanParametersExpanded'


class FixedParameterRange(ParameterRange):
    """Generic fixed parameter class.

    """

    def __init__(self, name, pRanged, pMin, pMax):
        super(FixedParameterRange, self).__init__(name, fixed=True)

        self.pRanged = pRanged
        self.pMin = pMin
        self.pMax = pMax
        self.name += '_%s_%f-%f' % (pRanged, pMin, pMax)


class DirectFixedParameters(FixedParameterRange):

    def __init__(self, name, pRanged, pIndex, pMin, pMax):
        super(DirectFixedParameters, self).__init__(name, pRanged, pMin, pMax)

        self.pIndex = pIndex
        self.pName = '%s_%d' % (pRanged, pIndex)
        self.pVal = 0.0
        self.name += '_curr_%f' % self.pVal
        self.direct = True


class SpecificParameterRange(ParameterRange):

    def __init__(self, name, pSet):
        super(SpecificParameterRange, self).__init__(name, specific=True)

        self.pSet = pSet


###############################################################################
# Classes for the paper
###############################################################################

class Figure1Parameters(FixedParameterRange):

    def __init__(self):
        super(Figure1Parameters, self).__init__('Figure1Parameters', 'pS',
              1.0e-06, 1.0)

        self.params = {'pR': (2.4e-02,), 'pS': (2.4e-03,),
                       'dR': (1e-04,), 'dS': (2.5e-04,),
                       'b': (1e-03,), 'u': (1e-03,),
                       'c': (0.035,), 'a': (0.5,)}

#class Figure1Parameters(SpecificParameterRange):

    #def __init__(self):
        #super(Figure1Parameters, self).__init__('Figure1Parameters', 'pS')

        #self.ps = [0.000024, 5.20652044412444E-05, 0.0001054032, 0.0004565523,
               #0.0009973093, 0.0019369488, 0.0024, 0.0026257313,
               #0.0026807684, 0.0029126971, 0.003208761, 0.0035105601,
               #0.0078836577, 0.0101821354, 0.0106134576, 0.0116117032,
               #0.0153114403, 0.018200827, 0.0236703803]

        #self.params = {'pR': (2.4e-02,), 'pS': (2.4e-2,),
                       #'dR': (1e-04,), 'dS': (2.5e-04,),
                       #'b': (1e-03,), 'u': (1e-03,),
                       #'c': (0.035,), 'a': (0.5,)}

        #self.name += '_pS'
#