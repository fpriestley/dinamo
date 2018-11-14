# DINAMO
DINAMO (DINAMO Is Not A MOdified blackbody) is a dust heating and emission code for astrophysics, 
described in Priestley, Barlow & De Looze (in prep).
Given dust optical properties, the local radiation field and the gas density and temperature, 
DINAMO calculates the equilibrium temperature distribution for each grain size in a specified 
size distribution, accounting for collisional heating by electrons and atoms/ions as well as 
radiative heating, stochastic heating effects and non-continuous cooling, 
and the resulting emitted spectral energy distribution. DINAMO has been benchmarked 
against other dust emission codes as described in Camps et al. (2015).

# User guide
1. Change the library path in the Makefile to point to an installation of the BLAS and LAPACK libraries.

2. Before compiling, check the hardcoded parameters are appropriate. The file constants_mod.f90 contains logical switches to change the behaviour of the code. Other potentially useful parameters are the number of temperature/enthalpy bins for each grain size (nEnthalpy in main.f90) and the minimum probability cutoff and the minimum temperature difference for convergence (minP and minTdiff in solve.f90).

3. Compile and run. The code reads input.in for the radiation and gas properties, and one file per grain species (specified in input.in) for the dust properties. Optical constants can be taken from the [Jena database](https://www.astro.uni-jena.de/Laboratory/OCDB/). To run the benchmarks, download the input data from [here](http://www.shg.ugent.be/html/index.html) and compile the code with lgBenchmark = true - the code will automatically read in the correct files from the input folder.
