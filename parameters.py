# -*- coding: utf-8 -*-
###############################################################################
# Classes for parameter ranges
###############################################################################
from __future__ import print_function


class ParameterRange(object):
    def __init__(self, name):
        self.name = name
        self.params = {}

    def get(self, param):
        if param in self.params:
            return self.params[param]
        else:
            print('Parameter %s not found!' % param)
            return None

    def set(self, param, val):
        self.params[param] = val


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
        self.set('pS', (2.4e-04, 2.4e-02))
        self.name = 'NitzanParametersReduced'


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