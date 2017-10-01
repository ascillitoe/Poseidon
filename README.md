# Poseidon

Poseidon - A high-order multi-block structured CFD code for executation on
           clusters of GPU's, many-core CPU's and accelerators.

Author(s) - Ashley Scillitoe

This is a multi-block structured CFD code using 6th order finite
differences and low storage RK to solve the compressible Navier-Stokes 
equations. The code is made parallel using the open source OPS 
framework (http://www.oerc.ox.ac.uk/projects/ops).

History:
16/09/17  - First build.                                             - as2341

I'm currently setting compiler and lib locations in libs/OPS/scripts/source_intellaptop. You'll need to change the absolute directories here. Probably a better way to do this but it is working for now!
