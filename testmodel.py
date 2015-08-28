# -*- coding: utf-8 -*-

import bngl

randParams = {}
randParams['vol'] = 2.0e-12              # cell volume (currently unused)
randParams['pR'] = (2.4e-03, 2.4e-01)    # transcription of R (kR)
randParams['pS'] = (2.4e-03, 2.4e-01)    # transcription of S (kS)
randParams['dR'] = (1e-05, 1e-03)        # decay of R (gR)
randParams['dS'] = (2.5e-05, 2.5e-03)    # decay of S (gS)
randParams['b'] = (1e-04, 1e-02)         # binding (k+)
randParams['u'] = (1e-04, 1e-02)         # unbinding (k-)
randParams['c'] = (7e-03, 7e-02)         # decay of complex (g)
randParams['a'] = (0.5, 0.5)             # alpha

m = bngl.CernetModel('/home/matt/testmodel', 1, 2, 1, 1, randParams)
