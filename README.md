# Poseidon

Poseidon - A high-order multi-block structured CFD code for executation on
           clusters of GPU's, many-core CPU's and accelerators.

Author(s) - Ashley Scillitoe

This is a multi-block structured CFD code using high-order finite
differences and low storage RK to solve the compressible Navier-Stokes 
equations. The code is made parallel using the open source OPS 
framework (http://www.oerc.ox.ac.uk/projects/ops).

History:
--------
16/09/17  - First build.                                             - as2341  
01/10/17  - Viscous terms added.                                     - as2341  


Build Instructions
------------------
1) OPS must first be installed (this needs the hdf5 library, MPI and CUDA libaries to be present, see https://github.com/gihanmudalige/OPS for more details).

2) Set compiler and lib enviroment variables. See src/sourcefiles/ for examples.

3) ./build.sh will run the OPS translator python scripts, and then compile the resulting code. 

(To compile in debug mode set DEBUG=true in build.sh)
