# bngl stuff for the finally final SimSys version

from __future__ import print_function, division
import os
import random
import networkx as nx
import matplotlib.pyplot as plt
import pexpect
from utilities import *


class CernetModel:
    def __init__(self, home, m, n, k, index, paramDict, logger, debugger,
                 seed=None, volScaling=False, template=None):

        if seed is not None:
            random.seed(seed)

        self.logger = logger
        self.debugger = debugger

        self.home = home
        self.plotDir = join(home, 'plots')
        self.name = 'ceRNET_%dx%d_%d_%d' % (m, n, k, index)
        #self.logger.info('Initializing model %d' % index)
        self.index = index
        self.filePath = join(home, self.name)
        self.bnglFile = self.filePath + '.bngl'
        makeDirs(self.plotDir)

        self.m = m
        self.n = n
        self.k = k
        self.ceRNAs = []
        self.miRNAs = []
        self.complexes = []
        self.corrMols = {}

        self.params = []
        self.molTypes = []
        self.observables = []
        self.rules = []

        if template is not None:
            self.ceRNAs = list(template.ceRNAs)
            self.miRNAs = list(template.miRNAs)
            self.molTypes = list(template.molTypes)
            self.complexes = list(template.complexes)
            self.observables = list(template.observables)

        self.createMolTypes()
        self.createComplexes()
        self.createObservables()
        self.createRulesAndParams(paramDict)

        self.writeNetworkFiles()

    def writeNetworkFiles(self):
        self.writeBNGL()
        self.writeSIF()
        self.writeGML()
        self.writeParamFile()
        self.saveNetworkImage()

    def copy(self, source, home, index, paramDict, seed):
        m = source.m
        n = source.n
        k = source.k
        logger = source.logger
        debugger = source.debugger
        c = CernetModel(home, m, n, k, index, paramDict, logger, debugger,
                        seed=seed, template=source)

        return c

    def getMolWithName(self, name):
        for mol in self.molTypes:
            if mol.name == name:
                return mol

        msg = 'CernetModel.findMolWithName did not find %s' % name
        self.debugger.warning(msg)
        return None

    def createMolTypes(self):
        count = 1
        for i in range(1, self.n + 1):
            self.ceRNAs.append(CeRNA(i, count))
            count += 1
        for i in range(1, self.m + 1):
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
            self.observables.append(
                BnglObservable('Molecules', '%s_free' % mol.name, mol.bnglCode))

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
            param = 'pT_%d' % mol.num
            self.params.append(BnglParameter(param, randVal))
            mol.prodRate = randVal

            # Create rule
            self.rules.append(BnglProductionRule(mol, param))

        for mol in self.miRNAs:
            # Choose param val
            minVal = paramDict['pR'][0]
            maxVal = paramDict['pR'][1]
            randVal = random.uniform(minVal, maxVal)
            param = 'pR_%d' % mol.num
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
            param = 'dT_%d' % mol.num
            self.params.append(BnglParameter(param, randVal))
            mol.decayRate = randVal

            # Create rule
            self.rules.append(BnglDecayRule(mol, param))

        for mol in self.miRNAs:
            # Choose param val
            minVal = paramDict['dR'][0]
            maxVal = paramDict['dR'][1]
            randVal = random.uniform(minVal, maxVal)
            param = 'dR_%d' % mol.num
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
            bName = 'b_{0}_{1}'.format(molX.num, molY.num)
            uName = 'u_{0}_{1}'.format(molX.num, molY.num)
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
            aName = 'a_{0}_{1}'.format(molX.num, molY.num)
            cName = 'c_{0}_{1}'.format(molX.num, molY.num)
            cFname = 'cF_{0}_{1}'.format(molX.num, molY.num)
            cPname = 'cP_{0}_{1}'.format(molX.num, molY.num)
            self.params.append(BnglParameter(aName, aRand))
            self.params.append(BnglParameter(cName, cRand))
            self.params.append(BnglParameter(cFname, '%s*%s' % (aName, cName)))
            self.params.append(BnglParameter(cPname,
                               '%s*(1-%s)' % (aName, cName)))

            miRNA = None
            if 'miRNA' in molX.name:
                miRNA = molX
            else:
                miRNA = molY

            # Full decay
            self.rules.append(BnglDecayRule(comp, cFname))
            # Partial decay
            self.rules.append(BnglDecayRule(comp, cPname,
                              remainder=miRNA.bnglCode))

    def writeBNGL(self):
        with open('%s.bngl' % self.filePath, 'w') as bnglFile:
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
        with open('%s.gml' % self.filePath, 'w') as gmlFile:
            gmlFile.write('graph [\n')
            gmlFile.write('\tlabel\t"%s"\n' % self.name)
            gmlFile.write('\tdirected\t0\n')

            for mol in self.ceRNAs:
                gmlFile.write('\tnode [\n')
                gmlFile.write('\t\tid\t%d\n' % mol.id)
                gmlFile.write('\t\tlabel\t"%s"\n' % mol.name)
                gmlFile.write('\t\ttype\t"%s"\n' % 'ceRNA')
                gmlFile.write('\t\tprodRate\t%f\n' % mol.prodRate)
                gmlFile.write('\t\tdecayRate\t%f\n' % mol.decayRate)
                gmlFile.write('\t\tgraphics [\n')
                gmlFile.write('\t\t\tfill\t"#0006ff"\n')
                gmlFile.write('\t\t]\t')
                gmlFile.write('\t]\n')

            for mol in self.miRNAs:
                gmlFile.write('\tnode [\n')
                gmlFile.write('\t\tid\t%d\n' % mol.id)
                gmlFile.write('\t\tlabel\t"%s"\n' % mol.name)
                gmlFile.write('\t\ttype\t"%s"\n' % 'miRNA')
                gmlFile.write('\t\tprodRate\t%f\n' % mol.prodRate)
                gmlFile.write('\t\tdecayRate\t%f\n' % mol.decayRate)
                gmlFile.write('\t\tgraphics [\n')
                gmlFile.write('\t\t\tfill\t"#ff0000"\n')
                gmlFile.write('\t\t]\t')
                gmlFile.write('\t]\n')

            for comp in self.complexes:
                gmlFile.write('\tedge [\n')
                gmlFile.write('\t\tsource\t%d\n' % comp.molX.id)
                gmlFile.write('\t\ttarget\t%d\n' % comp.molY.id)
                gmlFile.write('\t\tlabel\t"pp"\n')
                gmlFile.write('\t]\n')

            gmlFile.write(']\n')

    def writeSIF(self):
        with open('%s.sif' % self.filePath, 'w') as sifFile:
            for comp in self.complexes:
                s = '%s binds %s\n' % (comp.molX.name, comp.molY.name)
                sifFile.write(s)

    def writeParamFile(self):
        with open('%s.csv' % self.filePath, 'w') as csvFile:
            csvFile.write('parameter;value\n')
            for param in self.params:
                csvFile.write('%s;%s\n' % (param.name, param.val))

    def saveNetworkImage(self):
        gmlFile = self.filePath + '.gml'
        g = nx.read_gml(gmlFile, relabel=True)
        p = nx.circular_layout(g)
        nx.draw_networkx(g, pos=p)
        plt.savefig(self.filePath + '.png')
        plt.close()


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
        ts = '{0}({1}!1).{2}({3}!1)'
        self.bnglCode = ts.format(x.name, x.comps, y.name, y.comps)
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
        ts = '{0}({1}) + {2}({3}) <-> {0}({1}!1).{2}({3}!1) {4}, {5}'
        self.bnglCode = ts.format(x.name, x.comps, y.name, y.comps, b, u)


class BnglProductionRule(BnglRule):
    def __init__(self, x, k):
        self.mol = x
        self.k = k
        self.bnglCode = '0 -> {0} {1}'.format(x.bnglCode, k)


class BnglDecayRule(BnglRule):
    def __init__(self, x, k, remainder='0'):
        self.mol = x
        self.k = k
        ts = '{0} -> {1} {2} DeleteMolecules'
        self.bnglCode = ts.format(x.bnglCode, remainder, k)


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
    def __init__(self, model, logging=True, bngex=None):
        self.logging = logging
        if bngex is None:
            bngex = join(os.path.expanduser('~'), 'apps/BioNetGen/BNG2.pl')
        self.cmd = 'perl %s --console' % bngex
        self.model = model

        return

    def expect(self):
        self.console.expect('BNG>')
        return

    def loadModel(self, model):
        self.console.sendline('load ' + model)
        self.expect()
        return

    def initialize(self):
        self.console = pexpect.spawn(self.cmd, timeout=300000)
        if self.logging:
            fn = join(self.model.home, 'sim_log.log')
            self.console.logfile_read = open(fn, 'w')
        self.expect()
        self.loadModel(self.model.bnglFile)
        self.generateNetwork()
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
        self.console.close()
        return
