#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 09:19:49 2020

@author: thomas
"""

import numpy as np
from scipy.stats import norm

def Q(x):
    return norm.sf(x)

def symErrorQAM( SNRs, M ):    
    return 1.0 - (1.0 - 2.0 * Q( np.sqrt( 3.0 * SNRs / (M-1) ) ) ) ** 2.0

def bitErrorQAM( SNRs, M ):
    
    if M <= 1:
        return 1
    
    k = np.log2(M)
    return 1.0 / k * symErrorQAM(SNRs, M)

def makeArray(arr, num):
    if not isinstance(arr, np.ndarray):
        return np.ones( num ) * arr
    else:
        return arr
    
def subcarrierBER( eavs, N0s, Ms ):
    
    Pe = np.zeros( eavs.size )
    
    Ms = makeArray(Ms, eavs.size )
    N0s = makeArray(N0s, eavs.size )    
        
    for i, eav in enumerate(eavs):
        M = Ms[i]

        eav = eavs[i]
        N0 = N0s[i]
        SNR = eav / N0
        
        Pe[i] = bitErrorQAM ( SNR, M )
    
    return Pe

def averageSubcarrierBER( eavs, N0s, Ms ):
    totbits = np.sum( np.log2(Ms) )
    return np.sum( np.log2(Ms) * subcarrierBER(eavs, N0s, Ms) ) / totbits

def FoMAverageBER( eavs, fs, loss, N0s, ks, Pethr ):
    
    Ms = 2.0 ** ks
    eavs2 = eavs * loss(fs)
    Df = fs[1] - fs[0]    
    FoM = Df * np.sum( ks )    
    PeAvg = averageSubcarrierBER( eavs2, N0s, Ms )
    
    if PeAvg > Pethr:        
        return 0 
    else:
        return FoM

def FoMIndivdualBER( eavs, fs, loss, N0s, ks, Pethr):
    
    Ms = 2.0 ** ks
    eavs2 = eavs * loss(fs)
    Pe = subcarrierBER( eavs2, N0s, Ms)
    Df = fs[1] - fs[0]    
    
    return Df * np.sum( (Pe < Pethr)  * ks )
    
def FoM( params, fs, N0s, loss, Pethr, FoMtype = 'average'):
    N = int(params.size / 2)
    ks = params[0:N]
    eavs = params[N:]

    if FoMtype == 'average':        
        return FoMAverageBER( eavs, fs, loss, N0s, ks, Pethr )
    elif FoMtype == 'individual':        
        return FoMIndivdualBER( eavs, fs, loss, N0s, ks, Pethr )
    
def FoMflat( params, eavs, fs, N0s, loss, Pethr, FoMtype = 'average'):
    ks = params
    if FoMtype == 'average':        
        return FoMAverageBER( eavs, fs, loss, N0s, ks, Pethr )
    elif FoMtype == 'individual':        
        return FoMIndivdualBER( eavs, fs, loss, N0s, ks, Pethr )
    