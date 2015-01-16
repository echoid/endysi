# bngl stuff for the finally final SimSys version

import os
import sys
import math
import time
import glob
import random
from itertools import combinations
import networkx as nx
import pexpect
from utilities import *

class CernetModel:
    def __init__(self, name, rootDir, m, n, k, paramDict, seed=None, volScaling=False):
        if seed != None:
            random.seed(seed)
        
        self.home = home
        self.name = name
        self.path = join(home, name)
        makeDirs(self.home)
        self.m = m
        self.n = n
        self.k = k
        self.ceRNAs = []
        self.miRNAs = []
        self.complexes = []
        
        self.params = []
        self.molTypes = []
        self.observables = []
        self.rules = []
        
        self.createMolTypes()
        self.createComplexes()
        self.createObservables()
        self.createRulesAndParams(paramDict)
        
        self.writeBNGL()
        self.writeSIF()
        self.writeParamFile()
    
    
    def createMolTypes(self):
        count = 1
        for i in range(1, self.n+1):
            self.ceRNAs.append(CeRNA(i, count))
            count += 1
        for i in range(1, self.m+1):
            self.miRNAs.append(MiRNA(i, count))
            count += 1
        
        self.molTypes = list(self.ceRNAs)
        self.molTypes.extend(self.miRNAs)
    
    def createComplexes(self):
        for ceRNA in self.ceRNAs:
            regs = random.sample(self.miRNAs, self.k)
            ceRNA.regs.extend(regs)
            for reg in regs:
                reg.targs.append(ceRNA)
                self.complexes.append(BnglComplex(ceRNA, reg))
    
    def createObservables(self):
        for mol in self.molTypes:
            self.observables.append(BnglObservable('Molecules', '%s_free' % mol.name, mol.bnglCode))
    
    def createRulesAndParams(self, paramDict, volScaling=False):
        # Add basic cell params
        self.params.append(BnglParameter('V', paramDict['vol']))
        self.params.append(BnglParameter('NA', '6.0221415e+23'))
        self.params.append(BnglParameter('volScale', 'NA*V*1e-6'))
        
        # Production rules
        for mol in self.ceRNAs:
            # Choose param val
            minVal = paramDict['pT'][0]
            maxVal = paramDict['pT'][1]
            randVal = random.uniform(minVal, maxVal)
            param = 'pT%d' % mol.num
            self.params.append(BnglParameter(param, randVal))
            mol.prodRate = randVal
            
            # Create rule
            self.rules.append(BnglProductionRule(mol, param))
        
        for mol in self.miRNAs:
            # Choose param val
            minVal = paramDict['pR'][0]
            maxVal = paramDict['pR'][1]
            randVal = random.uniform(minVal, maxVal)
            param = 'pR%d' % mol.num
            self.params.append(BnglParameter(param, randVal))
            mol.prodRate = randVal
            
            # Create rule
            self.rules.append(BnglProductionRule(mol, param))
        
        # Decay rules
        for mol in self.ceRNAs:
            # Choose param val
            minVal = paramDict['dT'][0]
            maxVal = paramDict['dT'][1]
            randVal = random.uniform(minVal, maxVal)
            param = 'dT%d' % mol.num
            self.params.append(BnglParameter(param, randVal))
            mol.decayRate = randVal
            
            # Create rule
            self.rules.append(BnglDecayRule(mol, param))
        
        for mol in self.miRNAs:
            # Choose param val
            minVal = paramDict['dR'][0]
            maxVal = paramDict['dR'][1]
            randVal = random.uniform(minVal, maxVal)
            param = 'dR%d' % mol.num
            self.params.append(BnglParameter(param, randVal))
            mol.decayRate = randVal
            
            # Create rule
            self.rules.append(BnglDecayRule(mol, param))
        
        # Binding rules
        for comp in self.complexes:
            molX = comp.molX
            molY = comp.molY
            
            bMin = paramDict['b'][0]
            bMax = paramDict['b'][1]
            uMin = paramDict['u'][0]
            uMax = paramDict['u'][1]
        
            bRand = random.uniform(bMin, bMax)
            uRand = random.uniform(uMin, uMax)
            bName = 'b{0}{1}'.format(molX.num, molY.num)
            uName = 'u{0}{1}'.format(molX.num, molY.num)
            if volScaling:
                self.params.append(BnglParameter(bName, '%s/volScale' % bRand))
            else:
                self.params.append(BnglParameter(bName, bRand))
            self.params.append(BnglParameter(uName, uRand))
            comp.kon = bRand
            comp.koff = uRand
            
            # Create rule
            self.rules.append(BnglBindingRule(molX, molY, bName, uName))
        
        # Complex decay rules
        for comp in self.complexes:
            molX = comp.molX
            molY = comp.molY
            
            aMin = paramDict['a'][0]
            aMax = paramDict['a'][1]
            cMin = paramDict['c'][0]
            cMax = paramDict['c'][1]
            
            aRand = random.uniform(aMin, aMax)
            cRand = random.uniform(cMin, cMax)
            aName = 'a{0}{1}'.format(molX.num, molY.num)
            cName = 'c{0}{1}'.format(molX.num, molY.num)
            cFname = 'cF{0}{1}'.format(molX.num, molY.num)
            cPname = 'cP{0}{1}'.format(molX.num, molY.num)
            self.params.append(BnglParameter(aName, aRand))
            self.params.append(BnglParameter(cName, cRand))
            self.params.append(BnglParameter(cFname, '%s*%s' % (aName, cName)))
            self.params.append(BnglParameter(cPname, '%s*(1-%s)' % (aName, cName)))
            
            miRNA = None
            if 'miRNA' in molX.name:
                miRNA = molX
            else:
                miRNA = molY
            
            # Full decay
            self.rules.append(BnglDecayRule(comp, cFname))
            # Partial decay
            self.rules.append(BnglDecayRule(comp, cPname, remainder=miRNA.bnglCode))
    
    def writeBNGL(self):
        with open('%s.bngl' % self.path, 'w') as bnglFile:
            bnglFile.write('begin model\n')
            bnglFile.write('begin parameters\n')
            for param in self.params:
                bnglFile.write('\t%s\n' % param.bnglCode)
            bnglFile.write('end parameters\n')
            
            bnglFile.write('begin molecule types\n')
            for mol in self.molTypes:
                bnglFile.write('\t%s\n' % mol.bnglCode)
            bnglFile.write('end molecule types\n')
            
            bnglFile.write('begin observables\n')
            for obs in self.observables:
                bnglFile.write('\t%s\n' % obs.bnglCode)
            bnglFile.write('end observables\n')
            
            bnglFile.write('begin reaction rules\n')
            for rule in self.rules:
                bnglFile.write('\t%s\n' % rule.bnglCode)
            bnglFile.write('end reaction rules\n')
            
            bnglFile.write('end model\n')
    
    def writeGML(self):
        with open('%s.gml' % self.path, 'w') as gmlFile:
            gmlFile.write('Creator "Y"\n')
            gmlFile.write('Version 1.0\n')
            
            gmlFile.write('graph\n')
            gmlFile.write('[\n')
            gmlFile.write('\tlabel\t"%s"\n' % self.name)
            gmlFile.write('\tdirected\t0\n')
            
            for mol in self.ceRNAs:
                gmlFile.write('\tnode\n')
                gmlFile.write('\t\tid\t%d\n' % mol.id)
                gmlFile.write('\t\tlabel\t%s\n' % mol.name)
    
    def writeSIF(self):
        with open('%s.sif' % self.path, 'w') as sifFile:
            for comp in self.complexes:
                sifFile.write('%s binds %s\n' % (comp.molX.name, comp.molY.name))
    
    def writeParamFile(self):
        pass
    

class CeRNA:
    def __init__(self, num, id_):
        self.id = id_
        self.num = num
        self.name = 'ceRNA%d' % num
        self.comps = 'm'
        self.bnglCode = '%s(%s)' % (self.name, self.comps)
        self.regs = []
        self.prodRate = 0.0
        self.decayRate = 0.0

class MiRNA:
    def __init__(self, num, id_):
        self.id = id_
        self.num = num
        self.name = 'miRNA%d' % num
        self.comps = 'c'
        self.bnglCode = '%s(%s)' % (self.name, self.comps)
        self.targs = []
        self.prodRate = 0.0
        self.decayRate = 0.0

class BnglComplex:
    def __init__(self, x, y):
        self.molX = x
        self.molY = y
        self.bnglCode = '{0}({1}!1).{2}({3}!1)'.format(x.name, x.comps, y.name, y.comps)
        self.kon = 0.0
        self.koff = 0.0


class BnglRule:
    def __init__(self):
        pass
    

class BnglBindingRule(BnglRule):
    def __init__(self, x, y, b, u):
        self.molX = x
        self.molY = y
        self.kon = b
        self.koff = u
        self.bnglCode = '{0}({1}) + {2}({3}) <-> {0}({1}!1).{2}({3}!1) {4}, {5}'.format(x.name, x.comps, y.name, y.comps, b, u)
    

class BnglProductionRule(BnglRule):
    def __init__(self, x, k):
        self.mol = x
        self.k = k
        self.bnglCode = '0 -> {0} {1}'.format(x.bnglCode, k)


class BnglDecayRule(BnglRule):
    def __init__(self, x, k, remainder='0'):
        self.mol = x
        self.k = k
        self.bnglCode = '{0} -> {1} {2} DeleteMolecules'.format(x.bnglCode, remainder, k)
    
    
class BnglParameter:
    def __init__(self, param, val):
        self.name = param
        self.val = val
        self.bnglCode = '%s %s' % (param, val)

class BnglObservable:
    def __init__(self, type_, name, pattern):
        self.type = type_
        self.name = name 
        self.pattern = pattern
        self.bnglCode = '{0} {1} {2}'.format(type_, name, pattern)


class BnglSimulator:
    def __init__(self, modelFile=None, logging=False, bngex=None):
        # Create the console
        #plex = join(os.path.expanduser('~'), '.plenv/shims/perl')
        if bngex == None:
            bngex = join(os.path.expanduser('~'), 'apps/BioNetGen/BNG2.pl')
        cmd = 'perl %s --console' % bngex
        self.console = pexpect.spawn(cmd, timeout=30000)

        # Handle logging
        if logging:
            self.console.logfile_read = open(join(os.path.expanduser('~'), 'sim_log_read.log'), 'w')
            #self.console.logfile_send = sys.stdout
        self.expect()

        # Load model and initialize
        if not modelFile == None:
            self.initialize(modelFile)
        
        return
    
    def expect(self):
        self.console.expect('BNG>')
        return
    
    def loadModel(self, model):
        self.console.sendline('load ' + model)
        self.expect()
        return

    def initialize(self, modelFile):
        self.loadModel(modelFile)
        self.generateNetwork()
        self.equilibrate()
        return

    def generateNetwork(self):
        self.sendAction('generate_network({overwrite=>1})')
        return

    def saveConcentrations(self, label='initConcs'):
        self.sendAction('saveConcentrations("%s")' % label)
        return

    def resetConcentrations(self, label='initConcs'):
        self.sendAction('resetConcentrations("%s")' % label)
        return

    def setConcentration(self, species, conc):
        self.sendAction('setConcentration("%s","%s")' % (species, str(conc)))
        return

    def saveParameters(self, label='initParams'):
        self.sendAction('saveParameters("%s")' % label)
        return

    def resetParameters(self, label='initParams'):
        self.sendAction('resetParameters("%s")' % label)
        return

    def setParameter(self, param, val):
        self.sendAction('setParameter("%s","%s")' % (param, str(val)))
        return
    
    def sendAction(self, action):
        self.console.sendline('action ' + action)
        self.expect()
        return
    
    def clear(self):
        self.console.sendline('clear')
        self.expect()
        return

    def done(self):
        self.console.sendline('done')

    def close(self):
        #self.log.close()
        self.console.close()
        return

## testing
if __name__ == '__main__':
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
    
    cernet10by10 = CernetModel('c10by10', '/home/matt', 10, 10, 10, paramDict, seed=1)
    
    
