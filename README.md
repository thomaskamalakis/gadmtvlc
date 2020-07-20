# gadmtvlc

Genetic algorithm (GA)  optimization of a visible light communication using discrete multitone modulation
Implemented on Python 3.7

This project contains two sub-folders

analyticalmodel/
Provides an early version of the GA implementation where we used a preliminary VLC link model. File description
- mixed_ga.py : a library implementing a mixed integer genetic algorithm of Haupt, Randy L. "Antenna design with a mixed integer genetic algorithm." IEEE Transactions on Antennas and Propagation 55.3 (2007): 577-582.
- simpledmtlib.py : a library containing a simple model for the DMT/VLC system for evaluation
- dmtoptimization_mixedga.py : the "main" script demonstrating how the optimization can be carried out.

numericalmodel/
Provides an newer version of the GA implementation where we used a detailed VLC link model. File description
- mixed_ga.py : a library implementing a mixed integer genetic algorithm of Haupt, Randy L. "Antenna design with a mixed integer genetic algorithm." IEEE Transactions on Antennas and Propagation 55.3 (2007): 577-582. This is an improved version on what is found in the analyticalmodel/ folder containing logging and saving intermediate results
- libfom.py : a library containing some necessary functions for handling numerical fitness functions.
- libdmt.py : a lirbary containing the full blown VLC link model.
- runopt.py : the "main" script demonstrating how the optimization can be carried out.


