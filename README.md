# DINAMO
DINAMO (DINAMO Is Not A MOdified blackbody) is a dust heating and emission code for astrophysics, 
described in Priestley, Barlow & De Looze (2019).
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

# Input and output
Most of the input parameters and logical switches are fairly self-explanatory. The radiation field type must correspond to one of the functions defined in radfields.f90 - current options are blackbody, mathis (Mathis et al. 1983 ISM), pl (power law), sync (fit to Cassiopeia A synchrotron radiation) and crab (fit to Crab Nebula pulsar wind nebula emission). The dust optical data should be in the format from the Jena database: wavelength in microns or wavenumber in cm-1 (column 1), and real and imaginary parts of the refractive index (columns 2 and 3). Files using wavelength should end in '.lnk', files using wavenumber in '.vnk'. Although the files downloadable from the database are already in the correct format, **most of them do not extend to UV wavelengths so heating by the radiation field will not be properly calculated**. If radiative heating is going to be important, these data need to be extended to shorter wavelengths using e.g. the Draine & Lee (1984) data. Additionally, due to the wavelength grid being defined as that of the optical constants, **using multiple species will make bad things happen**. This should be fixed in a later version, but for now it's probably best to run one grain species at a time. The benchmark case works differently and runs all three species at once with no issues.

The code creates three files in output/, beginning with the model name. param.out contains the input parameters used for that model. The first line of Tdist.out is the number of grain species - after that, the file contains the species number and number of grain sizes, followed by the grain size in microns, temperature grid and probability distribution for each size. emis.out is the same as Tdist.out, except instead of temperature grid and probability it contains the emissivity in erg cm-3 s-1 sr-1 cm-1 for each grain size, followed by the wavelength grid in microns and the total emissivity from all grain sizes again in erg cm-3 s-1 sr-1 cm-1.
