# bcs_lattice_solver
Bardeen-Cooper-Schrieffer theory for superconductivity in a lattice Hubbard model at equilibrium specified by a temperature.
Solves the self-consistent superconducting gap and Hartree potential for an on-site Hubbard model with given attractive interaction strength.

## Installation
1. Get Eigen from eigen.tuxfamily.org and install hdf5 e.g. through Anaconda.
2. Link properly Eigen and hdf5 in the Makefile
3. Run `make` in the folder

 ## Running
 1. Specify the lattice geometry in the geometry.cfg file, following the example given in the file
 2. Specify other parameters and data points in parameters.cfg file, following the example given there.
 3. Specify the self-consistent algorithm chain in the scf_params.cfg file
 4. Run the program to generate data. Output is in hdf5 format, where the single-particle operator expectation values are included.
 5. Modify and run the given Python script analyze_data.py to plot the observables of interest
