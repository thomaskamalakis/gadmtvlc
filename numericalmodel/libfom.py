#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 14:01:08 2020

@author: thomas
"""
import sys
sys.path.insert(0,'/home/thomas/python/merlin2/merlin')

import random
import libdmt as dmt
import numpy as np

def sanityFun(p):
    return True

# figure of merit function
def FoM(k = None, scales = None, Df = 1e6, tapsize = 10, 
        nocarriers = 16, samplespersymbol = 10, 
        nobits = 1e4, cliplevel = 20, 
        psd = None, trf = None, Peth = 1e-3, 
        returnPhy = False,
        polynl = [0, 1, 0, 0]):

    M = 2 ** k
    M = M.astype(int)
    bitsPerCarrier = np.log2(M[1:])
    bitsPerFrame = np.sum(bitsPerCarrier)
    noFrames = np.ceil(nobits / bitsPerFrame).astype(int)
    
    np.random.seed(0)
    random.seed(0)
    
    phy=dmt.dmtphy(noframes = noFrames,
                   samplespersymbol = samplespersymbol,
                   tapsize = tapsize,
                   nocarriers = nocarriers,
                   Df = Df,
                   M = M,
                   polynl = polynl,
                   cliplevel = cliplevel,
                   psd = psd,
                   trfcallable = trf,
                   scales = scales)
    
    phy.digsimulate()
    
    
    if phy.BER < Peth:
        fom = phy.datarate        
    else:
        fom = 0
    
    np.random.seed()
    random.seed()
        
    if returnPhy:
        return fom, phy
    else:
        return fom