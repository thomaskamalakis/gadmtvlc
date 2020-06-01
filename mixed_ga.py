#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 09:46:17 2020

@author: thomas
"""

# The chromosome class

import numpy as np
import matplotlib.pyplot as plt

class chromosome:

    # initialize a chromosome
    # -initialGenes are the initial values of genes (can be omitted)
    # -sze is the number of genes
    # -mapType is a string 'RRRIIIRIR...' where 'R' indicates that the variable in 
    #          question is real and 'I' that it is integer, 'B' is for binary
    # -mins is an array containing the minimum values for each gene 
    # -maxs is an array containing the maximum values for each gene
    
    def __init__(self, 
                 initialGenes = None,
                 sze = None, 
                 mapType = None,
                 mins = None,
                 maxs = None,
                 fitnessFun = None,
                 sanityFun = None,
                 sanitize = True):
        
        if initialGenes is None:
            self.genes = np.random.rand( sze )
        else:
            self.genes = initialGenes
            
        self.mapType = mapType
        self.mins = mins
        self.maxs = maxs
        self.fitnessFun = fitnessFun
        self.sanityFun = sanityFun
        if sanitize:
            self.sanitize()
        
    def calcFitness( self ):
        values = self.variableValues()
        self.fitness = self.fitnessFun( values )
        
    def size( self ):
        return self.genes.size
    
    def sanitize(self):
        if self.sanityFun is not None:
            return self.sanityFun(self)
        else: 
            return True
    
    # Get actual variable Values    
    def variableValues(self):
        
        values = np.zeros(self.genes.size)
        
        for i, gene in enumerate(self.genes):
            
            if self.mapType[i] == 'R':
                values[i] = (self.maxs[i] - self.mins[i] ) * self.genes[i] \
                            + self.mins[i]        
                            
            elif self.mapType[i] == 'I':
                values[i] = (self.maxs[i] - self.mins[i] + 1) * self.genes[i] \
                            + self.mins[i]
                values[i] = np.floor( values[i] )
            
            elif self.mathType[i] == 'B':
                values[i] = np.round( values[i] )
                
        return values
    
    def geneValues(self,values):
        
        genes = np.zeros(values.size)
        
        for i, value in enumerate(values):
            
            if self.mapType[i] == 'R':
                genes[i] = ( values[i] - self.mins[i] ) / (self.maxs[i] - self.mins[i] )
            elif self.mapType[i] == 'I':
                genes[i] = ( values[i] - self.mins[i] ) / (self.maxs[i] - self.mins[i] + 1) 
        
        return genes
            
    def setGenesFromValues(self,values):
        self.genes = self.geneValues(values)
                
            
    # perform a correction factor mutation in the chromosome
    def mutateCorrectionFactor(self, mutationFactor = 0.1):
        
        ar = mutationFactor
        br = 1 - 2.0 * np.random.rand( self.genes.size ) 
        self.genes = mutationFactor * ar * (br * self.genes) + self.genes
        self.genes = self.genes - np.floor(self.genes)
        
    # general mutation method including correction factor, etc    
    def mutate(self, mutationType = 'correction factor', 
                     mutationFactor = 0.1, 
                     calcFitness = True,
                     sanitize = False):
        
        if mutationType == 'correction factor':
            self.mutateCorrectionFactor(mutationFactor = mutationFactor)
       
        if calcFitness:
            self.calcFitness()
        
        if sanitize:
            self.sanitize()
            
    def __str__(self):
        if hasattr(self, 'fitness'):
            return str('fitness = %e' %self.fitness)
        else:
            return str('fitness = None')
    
    def __repr__(self):
        return str('Chromsome with fitness = %e' %self.fitness)
        
    def check(self):
        values = self.variableValues()
        for i, v in enumerate(values):
            if not( self.mins[i] <= v <=self.maxs[i] ):
                return False
        
        return True
            
# Uniform crossover
def crossOverUniform(parent1, parent2, calcFitness = False, sanitize = False):
    
    c = chromosome(sze = parent1.size(), 
                   mapType = parent1.mapType,
                   mins = parent1.mins,
                   maxs = parent1.maxs,
                   fitnessFun = parent1.fitnessFun,
                   sanityFun = parent1.sanityFun,
                   sanitize = False)
    
    bMap = np.random.randint(0, high = 2, size = c.size() )
    
    for i, b in enumerate(bMap):
        
        if b == 0:
            c.genes[i] = parent1.genes[i]
        elif b == 1:
            c.genes[i] = parent2.genes[i]
    
    if calcFitness:
        c.calcFitness()
    
    if sanitize:
        c.sanitize()
        
    return c
        

class population:
    
    def __init__(self, 
                 noChromosomes = 20,
                 noGenes = 10,
                 initPopulation = None,
                 mapType = None,
                 mins = None,
                 maxs = None,
                 fitnessFun = None,
                 maxNoCrossovers = 100,
                 mutationFactor = 0.1,
                 verboseLvl = 0,
                 sanityFun = None):
        
        self.mapType = mapType
        self.mins = mins
        self.maxs = maxs
        self.fitnessFun = fitnessFun
        self.chromosomes = []
        self.fitnesses = []
        self.maxNoCrossovers = maxNoCrossovers
        self.crossovers = 0
        self.fitnessEvaluations = 0
        self.mutationFactor = mutationFactor
        self.verboseLvl = verboseLvl
        self.sanityFun = sanityFun
        
        self.bestFitnessRecords = []
        self.worstFitnessRecords = []
        
        if initPopulation is not None:
            self.noChromosomes, self.noGenes = initPopulation.shape
            for i in range(0, noChromosomes):
                c = chromosome(initialGenes = initPopulation[i,:],
                           mapType = self.mapType,
                           mins = self.mins, maxs = self.maxs, 
                           fitnessFun = self.fitnessFun,
                           sanityFun = self.sanityFun)
                 
                c.calcFitness()
                self.fitnessEvaluations += 1
                self.chromosomes.append( c )
                
                
        else:
            self.noChromosomes = noChromosomes
            self.noGenes = noGenes
            for i in range(0, noChromosomes):
                c = chromosome(sze = noGenes,
                           mapType = self.mapType,
                           mins = self.mins, maxs = self.maxs, 
                           fitnessFun = self.fitnessFun,
                           sanityFun = self.sanityFun) 
                c.calcFitness()
                self.fitnessEvaluations += 1                
                self.chromosomes.append( c )
                
        self.sort()
    
    def getFitnesses(self):
        return np.array( [x.fitness for x in self.chromosomes] )
            
    def sort(self):
        self.chromosomes.sort(key = lambda x: x.fitness, reverse = True)
        
    def tournamentSelection(self, tournamentSize = 2, level = 0.5 ):
        
        pool = []
        sze = tournamentSize
        self.sort()
        poolSize = round(level * self.noChromosomes)
        
        while sze >= 1:
            r = np.random.randint(0, poolSize, 1)[0]
            if r not in pool:
                pool.append(r)
                sze -= 1
                
        return pool
    
    def select(self, selectionType = 'tournament', tournamentSize = 2, level = 0.5):
        
        if selectionType == 'tournament':
            return self.tournamentSelection(tournamentSize = tournamentSize)
            
    def crossOver(self, i, j, crossOverType = 'uniform', sanitize = False):
        if crossOverType =='uniform':
           return crossOverUniform( self.chromosomes[i], self.chromosomes[j], sanitize = sanitize)
        self.crossovers += 1
        
    def replace(self, chromosome):
        
        fitnesses = self.getFitnesses()
        if chromosome.fitness > np.min( fitnesses ):
            i = np.argmin( fitnesses )
            self.chromosomes[i] = chromosome
        
        self.sort()
            
    def __str__(self):
        st = 'population with %d chromosomes and %d genes per chromosome\n' %(self.noChromosomes, self.noGenes)
        st = st + 'number of crossovers: %d \n' %self.crossovers
        st = st + 'number of fitness evaluations: %d\n' %self.fitnessEvaluations
        st = st + 'best fitness so far: %e\n' %self.chromosomes[0].fitness

        return st
    
    def updateRecords(self):
        self.bestFitnessRecords.append( self.chromosomes[0].fitness )
        self.worstFitnessRecords.append( self.chromosomes[self.noChromosomes - 1].fitness )
        
        
    def simulate(self):
        self.crossovers = 0
        
        while (self.crossovers < self.maxNoCrossovers):
            if self.verboseLvl >= 1:
                print(self)
            pool = self.tournamentSelection()
            i = pool[0]
            j = pool[1]
            if self.verboseLvl >= 2:
                print('Selecting chromosomes with fitness %e and %e for pairing' %(self.chromosomes[i].fitness, self.chromosomes[j].fitness) )
            ch = self.crossOver(i, j)
            self.crossovers += 1                
            if self.verboseLvl >= 2:            
                print('Resulting chromosome')
            
            
            ch.mutate(mutationFactor = self.mutationFactor, sanitize = True)            
            self.fitnessEvaluations += 1
            
            if self.verboseLvl > 0:
                print(ch)            
            self.replace(ch)
            self.updateRecords()
            
