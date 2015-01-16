####################################################
# Pexpect-based interface for the BNG console.
#
####################################################


import pexpect
import sys

class BNGCI:
    def __init__(self, modelFile=None, logging=False):
        # Create the console
        self.console = pexpect.spawn('BNG2.pl -console', timeout=300000)

        # Handle logging
        if logging:
            self.console.logfile_read = file('sim_log_read.log', 'w')
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
        self.sendAction('setConcentration("%s",%s)' % (species, str(conc)))
        return

    def saveParameters(self, label='initParams'):
        self.sendAction('saveParameters("%s")' % label)
        return

    def resetParameters(self, label='initParams'):
        self.sendAction('resetParameters("%s")' % label)
        return

    def setParameter(self, param, val):
        self.sendAction('setParameter("%s",%s)' % (param, str(val)))
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
        self.sendAction('done')

    def close(self):
        #self.log.close()
        self.console.close()
        return

    ##### testing #####
    def __setattr__(self, name, value):
        self.__dict__[name] = value
        return
    
    

    
