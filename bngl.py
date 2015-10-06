# Core classes representing BioNetGen elements and ceRNA network components

from __future__ import print_function, division
import os
import random
import math
import networkx as nx
import matplotlib.pyplot as plt
import pexpect
from utilities import *

# Global params:
_nFlips = 10000

# Debugging
_debugging = False


class CernetModel:
    def __init__(self, home, m, n, k, index, paramRange, seed=None, alpha=None,
                 volScaling=False, template=None, linearSampling=False):

        if seed is not None:
            random.seed(seed)
        self.seed = seed

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
        self.alpha = alpha
        self.ceRNAs = []
        self.miRNAs = []
        self.complexes = []
        self.corrMols = {}

        self.params = []
        self.molTypes = []
        self.observables = []
        self.rules = []

        count = 1
        count = self.createMolTypes(count)
        count = self.createComplexes(count)
        self.createObservables()
        self.createRulesAndParams(paramRange, linearSampling=linearSampling,
                                  alpha=alpha)
        self.writeNetworkFiles()

    def purge(self):
        del self.ceRNAs
        del self.miRNAs
        del self.complexes
        del self.corrMols
        del self.params
        del self.molTypes
        del self.observables
        del self.rules

    def writeNetworkFiles(self):
        self.writeBNGL()
        self.writeNet()
        self.writeSIF()
        self.writeGML()
        self.writeParamFile()
        #self.saveNetworkImage()

    def getMolWithName(self, name):
        for mol in self.molTypes:
            if mol.name == name:
                return mol

        msg = 'CernetModel.findMolWithName did not find %s' % name
        self.debugger.warning(msg)
        return None

    def createMolTypes(self, count):
        for i in range(1, self.n + 1):
            self.ceRNAs.append(CeRNA(i, count))
            count += 1
        for i in range(1, self.m + 1):
            self.miRNAs.append(MiRNA(i, count))
            count += 1

        self.molTypes = list(self.ceRNAs)
        self.molTypes.extend(self.miRNAs)

        return count

    def createComplexes_orig(self, count):
        for ceRNA in self.ceRNAs:
            partners = random.sample(self.miRNAs, self.k)
            ceRNA.partners.extend(partners)
            for reg in partners:
                reg.partners.append(ceRNA)
                self.complexes.append(BnglComplex(ceRNA, reg, count))
                count += 1

        return count

    def createComplexes(self, count):
        if self.m <= 1 or self.k == 0 or self.n == self.k == self.m:
            self.createComplexes_orig(count)

        else:
            for ceRNA in self.ceRNAs:
                # get the indices for the miRNAs
                indices = []
                for i in range(self.k):
                    d = ceRNA.num + i
                    if d > self.m:
                        d = d - self.m
                    indices.append(d)

                # create the connections
                for i in indices:
                    miRNA = self.getMolWithName('miRNA%d' % i)
                    self.complexes.append(BnglComplex(ceRNA, miRNA, count))
                    ceRNA.partners.append(miRNA)
                    miRNA.partners.append(ceRNA)
                    count += 1

            self.shuffleNodes()

            return count

    def shuffleNodes(self):
        for i in range(_nFlips):
            comps = random.sample(self.complexes, 2)

            while not self.swapNodes(comps[0], comps[1]):
                comps = random.sample(self.complexes, 2)

    def swapNodes(self, compA, compB):
        ax = compA.molX
        ay = compA.molY
        bx = compB.molX
        by = compB.molY

        if ax.hasPartner(by) or bx.hasPartner(ay):
            return False

        ax.replacePartner(ay, by)
        ay.replacePartner(ax, bx)
        bx.replacePartner(by, ay)
        by.replacePartner(bx, ax)

        compA.molX = bx
        compB.molX = ax
        compA.rewriteCode()
        compB.rewriteCode()
        return True

    def createObservables(self):
        count = 1
        for mol in self.molTypes:
            self.observables.append(
                BnglObservable('Molecules', '%s_free' % mol.name,
                                mol, count))
            count += 1

    def createRulesAndParams(self, paramRange, volScaling=False, alpha=None,
                             linearSampling=False):
        # Add basic cell params
        #self.params.append(BnglParameter('V', paramRange['vol']))
        #self.params.append(BnglParameter('NA', '6.0221415e+23'))
        #self.params.append(BnglParameter('volScale', 'NA*V*1e-6'))

        count = 0
        pCount = rCount = 1
        pVal = 0.0
        # Production rules
        for mol in self.ceRNAs:
            # Choose param val
            if paramRange.isFixed():
                pVal = paramRange.getMin('pR')
            else:
                minVal = paramRange.getMin('pR')
                maxVal = paramRange.getMax('pR')
                if linearSampling:
                    pVal = random.uniform(minVal, maxVal)
                else:
                    pVal = math.exp(random.uniform(math.log(minVal),
                                                  math.log(maxVal)))
            param = 'pR_%d' % mol.num
            self.params.append(BnglParameter(param, pVal, pCount))
            mol.prodRate = pVal

            # Create rule
            self.rules.append(BnglProductionRule(mol, param, rCount))
            count += 1
            pCount += 1
            rCount += 1

        # confirm correct number of rules made
        assert count == self.n

        count = 0
        for mol in self.miRNAs:
            # Choose param val
            if paramRange.isFixed():
                pVal = paramRange.getMin('pS')
            else:
                minVal = paramRange.getMin('pS')
                maxVal = paramRange.getMax('pS')
                if linearSampling:
                    pVal = random.uniform(minVal, maxVal)
                else:
                    pVal = math.exp(random.uniform(math.log(minVal),
                                                  math.log(maxVal)))
            param = 'pS_%d' % mol.num
            self.params.append(BnglParameter(param, pVal, pCount))
            mol.prodRate = pVal

            # Create rule
            self.rules.append(BnglProductionRule(mol, param, rCount))
            count += 1
            pCount += 1
            rCount += 1

        assert count == self.m

        # Decay rules
        count = 0
        for mol in self.ceRNAs:
            # Choose param val
            if paramRange.isFixed():
                pVal = paramRange.getMin('dR')
            else:
                minVal = paramRange.getMin('dR')
                maxVal = paramRange.getMax('dR')
                if linearSampling:
                    pVal = random.uniform(minVal, maxVal)
                else:
                    pVal = math.exp(random.uniform(math.log(minVal),
                                                      math.log(maxVal)))
            param = 'dR_%d' % mol.num
            self.params.append(BnglParameter(param, pVal, pCount))
            mol.decayRate = pVal

            # Create rule
            self.rules.append(BnglDecayRule(mol, param, rCount))
            count += 1
            pCount += 1
            rCount += 1

        assert count == self.n

        count = 0
        for mol in self.miRNAs:
            # Choose param val
            if paramRange.isFixed():
                pVal = paramRange.getMin('dS')
            else:
                minVal = paramRange.getMin('dS')
                maxVal = paramRange.getMax('dS')
                if linearSampling:
                    pVal = random.uniform(minVal, maxVal)
                else:
                    pVal = math.exp(random.uniform(math.log(minVal),
                                                      math.log(maxVal)))
            param = 'dS_%d' % mol.num
            self.params.append(BnglParameter(param, pVal, pCount))
            mol.decayRate = pVal

            # Create rule
            self.rules.append(BnglDecayRule(mol, param, rCount))
            count += 1
            pCount += 1
            rCount += 1

        assert count == self.m

        count = 0
        for comp in self.complexes:
            molX = comp.molX
            molY = comp.molY

            if paramRange.isFixed():
                bRand = paramRange.getMin('b')
                uRand = paramRange.getMin('u')
            else:
                bMin = paramRange.getMin('b')
                bMax = paramRange.getMax('b')
                uMin = paramRange.getMin('u')
                uMax = paramRange.getMax('u')

                if linearSampling:
                    bRand = random.uniform(bMin, bMax)
                    uRand = random.uniform(uMin, uMax)
                else:
                    bRand = math.exp(random.uniform(math.log(bMin),
                                                    math.log(bMax)))
                    uRand = math.exp(random.uniform(math.log(uMin),
                                                    math.log(uMax)))

            bName = 'b_{0}_{1}'.format(molX.num, molY.num)
            uName = 'u_{0}_{1}'.format(molX.num, molY.num)

            if volScaling:
                self.params.append(BnglParameter(bName, '%s/volScale' % bRand,
                                   pCount))
            else:
                self.params.append(BnglParameter(bName, bRand, pCount))

            self.params.append(BnglParameter(uName, uRand, pCount + 1))
            comp.kon = bRand
            comp.koff = uRand

            # Create rule
            self.rules.append(BnglBindingRule(comp, bName, uName, rCount))
            count += 1
            pCount += 2
            rCount += 2

        assert count == (self.n * self.k)

        # Complex decay rules
        count = 0
        for comp in self.complexes:
            molX = comp.molX
            molY = comp.molY

            if paramRange.isFixed():
                aRand = paramRange.getMin('a')
                cRand = paramRange.getMin('c')
            else:
                aMin = paramRange.getMin('a')
                aMax = paramRange.getMax('a')
                cMin = paramRange.getMin('c')
                cMax = paramRange.getMax('c')

                if linearSampling:
                    if alpha is None:
                        aRand = random.uniform(aMin, aMax)
                    else:
                        aRand = alpha

                    cRand = random.uniform(cMin, cMax)
                else:
                    if alpha is None:
                        aRand = math.exp(random.uniform(math.log(aMin),
                                                        math.log(aMax)))
                    else:
                        aRand = alpha

                    cRand = math.exp(random.uniform(math.log(cMin),
                                                    math.log(cMax)))

            aName = 'a_{0}_{1}'.format(molX.num, molY.num)
            cName = 'c_{0}_{1}'.format(molX.num, molY.num)
            cFname = 'cF_{0}_{1}'.format(molX.num, molY.num)
            cPname = 'cP_{0}_{1}'.format(molX.num, molY.num)
            self.params.append(BnglParameter(aName, aRand, pCount))
            self.params.append(BnglParameter(cName, cRand, pCount + 1))
            self.params.append(BnglParameter(cFname, '%s*%s' % (aName, cName),
                                             pCount + 2))
            self.params.append(BnglParameter(cPname,
                               '%s*(1-%s)' % (aName, cName), pCount + 3))

            miRNA = None
            if 'miRNA' in molX.name:
                miRNA = molX
            else:
                miRNA = molY

            # Full decay
            self.rules.append(BnglDecayRule(comp, cFname, rCount))
            # Partial decay
            self.rules.append(BnglDecayRule(comp, cPname, rCount + 1,
                              remainder=miRNA))
            count += 1
            pCount += 4
            rCount += 2

        assert count == (self.n * self.k)

    def createRulesAndParams_log(self, paramRange, volScaling=False):
        # Production rules
        count = 0
        for mol in self.ceRNAs:
            # Choose param val
            minVal = paramRange['pR'][0]
            maxVal = paramRange['pR'][1]
            #pVal = 10.0 ** math.log10(random.uniform(minVal, maxVal))
            pVal = math.exp(random.uniform(math.log(minVal),
                                              math.log(maxVal)))
            param = 'pR_%d' % mol.num
            self.params.append(BnglParameter(param, pVal))
            mol.prodRate = pVal

            # Create rule
            self.rules.append(BnglProductionRule(mol, param))
            count += 1

        # confirm that the right number of rules are created
        assert count == self.n

        count = 0
        for mol in self.miRNAs:
            # Choose param val
            minVal = paramRange['pS'][0]
            maxVal = paramRange['pS'][1]
            pVal = math.exp(random.uniform(math.log(minVal),
                                              math.log(maxVal)))
            param = 'pS_%d' % mol.num
            self.params.append(BnglParameter(param, pVal))
            mol.prodRate = pVal

            # Create rule
            self.rules.append(BnglProductionRule(mol, param))
            count += 1

        assert count == self.m

        # Decay rules
        count = 0
        for mol in self.ceRNAs:
            # Choose param val
            minVal = paramRange['dR'][0]
            maxVal = paramRange['dR'][1]
            pVal = math.exp(random.uniform(math.log(minVal),
                                              math.log(maxVal)))

            param = 'dR_%d' % mol.num
            self.params.append(BnglParameter(param, pVal))
            mol.decayRate = pVal

            # Create rule
            self.rules.append(BnglDecayRule(mol, param))
            count += 1

        assert count == self.n

        count = 0
        for mol in self.miRNAs:
            # Choose param val
            minVal = paramRange['dS'][0]
            maxVal = paramRange['dS'][1]
            pVal = math.exp(random.uniform(math.log(minVal),
                                              math.log(maxVal)))

            param = 'dS_%d' % mol.num
            self.params.append(BnglParameter(param, pVal))
            mol.decayRate = pVal

            # Create rule
            self.rules.append(BnglDecayRule(mol, param))
            count += 1

        assert count == self.m

        count = 0
        for comp in self.complexes:
            molX = comp.molX
            molY = comp.molY

            bMin = paramRange['b'][0]
            bMax = paramRange['b'][1]
            uMin = paramRange['u'][0]
            uMax = paramRange['u'][1]

            bRand = math.exp(random.uniform(math.log(bMin),
                                              math.log(bMax)))
            uRand = math.exp(random.uniform(math.log(uMin),
                                              math.log(uMax)))

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
            count += 1

        assert count == (self.n * self.k)

        # Complex decay rules
        count = 0
        for comp in self.complexes:
            molX = comp.molX
            molY = comp.molY

            aMin = paramRange['a'][0]
            aMax = paramRange['a'][1]
            cMin = paramRange['c'][0]
            cMax = paramRange['c'][1]

            aRand = math.exp(random.uniform(math.log(aMin),
                                              math.log(aMax)))
            cRand = math.exp(random.uniform(math.log(cMin),
                                              math.log(cMax)))

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
                              remainder=miRNA))

            count += 1

        assert count == (self.n * self.k)

        return

    def writeNet(self):
        with open('%s.net' % self.filePath, 'w') as netFile:
            netFile.write('begin parameters\n')
            for param in self.params:
                netFile.write('\t%s\n' % param.netCode)
            netFile.write('end parameters\n')

            netFile.write('begin species\n')
            for mol in self.molTypes:
                netFile.write('\t%d %s 0\n' % (mol.id, mol.bnglCode))
            for comp in self.complexes:
                netFile.write('\t%d %s 0\n' % (comp.id, comp.bnglCode))
            netFile.write('end species\n')

            netFile.write('begin reactions\n')
            for rule in self.rules:
                netFile.write('\t%d %s\n' % (rule.id, rule.netCode))
            netFile.write('end reactions\n')

            netFile.write('begin groups\n')
            for obs in self.observables:
                netFile.write('\t%s\n' % obs.netCode)
            netFile.write('end groups')

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


class RNA:
    def __init__(self):
        self.partners = []

    def hasDuplicates(self, r):
        return r in self.partners

    def hasPartner(self, p):
        return p in self.partners

    def replacePartner(self, oldP, newP):
        if oldP not in self.partners:
            return

        self.partners.remove(oldP)
        self.partners.append(newP)


class CeRNA(RNA):
    def __init__(self, num, id_):
        self.id = id_
        self.num = num
        self.name = 'ceRNA%d' % num
        self.comps = 'm'
        self.bnglCode = '%s(%s)' % (self.name, self.comps)
        self.netCode = str(self.id)
        self.partners = []
        self.prodRate = 0.0
        self.decayRate = 0.0


class MiRNA(RNA):
    def __init__(self, num, id_):
        self.id = id_
        self.num = num
        self.name = 'miRNA%d' % num
        self.comps = 'c'
        self.bnglCode = '%s(%s)' % (self.name, self.comps)
        self.netCode = str(self.id)
        self.partners = []
        self.prodRate = 0.0
        self.decayRate = 0.0


class BnglComplex:
    def __init__(self, x, y, id_):
        self.id = id_
        self.molX = x
        self.molY = y
        ts = '{0}({1}!1).{2}({3}!1)'
        self.bnglCode = ts.format(x.name, x.comps, y.name, y.comps)
        self.netCode = str(self.id)
        self.kon = 0.0
        self.koff = 0.0
        self.decayRate = 0.0

    def rewriteCode(self):
        ts = '{0}({1}!1).{2}({3}!1)'
        self.bnglCode = ts.format(self.molX.name, self.molX.comps,
                                  self.molY.name, self.molY.comps)


class BnglRule:
    def __init__(self):
        self.reactants = []
        self.products = []
        pass

    def hasReactant(self, x):
        return x in self.reactants

    def hasProduct(self, x):
        return x in self.products


class BnglBindingRule(BnglRule):
    def __init__(self, comp, b, u, id_):
        self.id = id_
        x = comp.molX
        y = comp.molY
        self.comp = comp
        ts = '{0}({1}) + {2}({3}) <-> {0}({1}!1).{2}({3}!1) {4}, {5}'
        self.bnglCode = ts.format(x.name, x.comps, y.name, y.comps, b, u)
        ts = '{0},{1} {2} {3}\n\t{5} {2} {0},{1} {4}'
        self.netCode = ts.format(x.id, y.id, comp.id, b, u, id_ + 1)


class BnglProductionRule(BnglRule):
    def __init__(self, x, k, id_):
        self.id = id_
        self.mol = x
        self.k = k
        self.bnglCode = '0 -> {0} {1}'.format(x.bnglCode, k)
        self.netCode = '0 {0} {1}'.format(x.netCode, k)


class BnglDecayRule(BnglRule):
    def __init__(self, x, k, id_, remainder='0'):
        self.id = id_
        self.mol = x
        self.k = k
        self.remainder = remainder
        r = rn = '0'
        if remainder != '0':
            r = remainder.bnglCode
            rn = remainder.netCode
        ts = '{0} -> {1} {2} DeleteMolecules'
        self.bnglCode = ts.format(x.bnglCode, r, k)
        self.netCode = '{0} {1} {2}'.format(x.netCode, rn, k)


class BnglParameter:
    def __init__(self, param, val, id_):
        self.id = id_
        self.name = param
        self.val = val
        self.bnglCode = '%s %s' % (param, val)
        self.netCode = '{0} {1} {2}'.format(id_, param, val)


class BnglObservable:
    def __init__(self, type_, name, mol, id_):
        self.id = id_
        self.type = type_
        self.name = name
        self.mol = mol
        self.bnglCode = '{0} {1} {2}'.format(type_, name, mol.bnglCode)
        self.netCode = '{0} {1} {2}'.format(id_, mol.name, mol.netCode)


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
