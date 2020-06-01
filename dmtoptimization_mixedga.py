#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:38:31 2020

@author: thomas
"""

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
import mixed_ga as ga
from simpledmtlib import FoMflat
#
#def sanityFun(c):
#    values = c.variableValues()
#    N = int(values.size / 2)
#    powerDis = values[N:]
#    powerDis = powerDis / np.sum(powerDis)
#    values[N:] = powerDis
#    c.setGenesFromValues(values)    
#    return True     

def sanityFun(c):
    return True
# Subcarrier number    
N = 16

# Maximum carrier frequency allowed
Fmax = 1e6

# Channel pwoer transfer function
loss = lambda f: ( 1- 0.5 * f/Fmax  ) * 2

# BER threshold 
Peth = 1e-3


Df = 1e6 / N
fs = np.arange(0, N) * Df
Ts = 1 / Df
Peth = 1e-3
nbits = 3
eavs = loss(fs)

N0 = 1e-2
Pethr = 1e-3
plt.close('all')

print('Setting Up Population')
mapType = 'I' * N 
mins = 2 * np.ones(N)
maxs = 8 * np.ones(N)

np.random.seed(1)

fitnessFun = lambda params :  FoMflat( params, eavs, fs, N0, loss, Pethr, FoMtype = 'individual') / 1e6

p = ga.population(noChromosomes = 50,
                  noGenes = N,
                  mapType = mapType,
                  mins = mins,
                  maxs = maxs,
                  fitnessFun = fitnessFun,
                  maxNoCrossovers = 2400, 
                  mutationFactor = 0.1,
                  verboseLvl = 1,
                  sanityFun = sanityFun)
#        
p.simulate()
values = p.chromosomes[0].variableValues()
ks = values[0:N]
plt.rc('text', usetex=True)
matplotlib.rc('lines', linewidth=2)
matplotlib.rc('font', family='sans-serif') 
matplotlib.rc('font', serif='Helvetica') 
matplotlib.rcParams.update({'font.size': 14})
plt.close('all')
plt.figure()
plt.plot( p.bestFitnessRecords, label = "$\mathrm{max}\{R'_\mathrm{b}\}$" )
plt.plot( p.worstFitnessRecords,'--', label = "$\mathrm{min}\{R'_\mathrm{b}\}$" )
plt.yticks((2, 3, 4, 5, 6))
plt.ylabel("$R'_\mathrm{b}$")
plt.xlabel('$j$')
plt.legend()
plt.figure()
plt.plot(ks)
plt.yticks((2, 3, 4, 5, 6, 7))
plt.ylabel("$\log_2M_k$")
plt.xlabel('$k$')


