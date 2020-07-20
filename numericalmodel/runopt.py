#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:38:31 2020

@author: thomas
"""

import numpy as np
import matplotlib.pyplot as plt
import mixed_ga as ga

from libfom import FoM, sanityFun

#
#def sanityFun(c):
#    values = c.variableValues()
#    N = int(values.size / 2)
#    powerDis = values[N:]
#    powerDis = powerDis / np.sum(powerDis)
#    values[N:] = powerDis
#    c.setGenesFromValues(values)    
#    return True     

plt.close('all')
# Subcarrier number    
N = 16

# Maximum carrier frequency allowed
Fmax = 3.0e6

# Subcarrier spacing
Df = Fmax / N

# LED 3dB frequency
f0 = 1e6

# Channel pwoer transfer function
trf = lambda x: 1.0/(1j*x/f0+1)

# noise psd
Hf = lambda f: 0.5e-3

# BER threshold 
Peth = 1e-3

# number of variables
Nparams = 2 * N + 2 

# type of parameters
mapType = 'I' * N + 'R' * (N+2)

# bounds for parameters
mins = 2 * np.ones(Nparams)
maxs = 6 * np.ones(Nparams)

mins[N] = 8
maxs[N] = 15

mins[N+1] = 0.25e6
maxs[N+1] = Fmax

mins[N+2:] = 0.5
maxs[N+2:] = 1.0

nobits = 5e4
tapsize = 10   
polynl = [0, 1, 0.01, 0.02] 

   
fitnessFun = lambda params :  FoM( scales = params[N+2:], k = params[0:N].astype(int), 
                                   Df = params[N+1] / (N-1), 
                                   tapsize = tapsize, 
                                   nocarriers = N, 
                                   nobits = nobits, 
                                   cliplevel = params[N], 
                                   psd = Hf, 
                                   trf = trf,
                                   polynl = polynl)

p = ga.population(noChromosomes = 400,
                  noGenes = Nparams,
                  mapType = mapType,
                  mins = mins,
                  maxs = maxs,
                  fitnessFun = fitnessFun,
                  maxNoCrossovers = 20000, 
                  saveData = True,
                  saveEvery = 1000,
                  mutationFactor = 0.2,
                  verboseLvl = 2,
                  filename = 'optsctracked',                  
                  sanityFun = sanityFun,
                  keepVariableTrack = True,
                  seed = 0)
        
p.simulate()

plt.figure(1)
plt.plot( p.bestFitnessRecords )
plt.plot( p.worstFitnessRecords )

params = p.chromosomes[0].variableValues()
k = params[0:N]
cliplevel = params[N]
Fmax = params[N+1]
M = 2**k
M = M.astype(int)
scales = params[N+2:]

