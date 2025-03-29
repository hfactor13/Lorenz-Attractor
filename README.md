# About

The purpose of this repository is to implement a common numerical algorithm to solve differential equations using the Fortran language. The calculations from Fortran are then compared with the Python implementation.

# Getting Started

1. Execute the following command on the command line `gfortran lorenz_attractor.f90 linspace_module.f90 -o ./output/lorenz_attractor` (make sure to create the output directory in advance)
2. Execute the executable `./output/lorenz_attractor` to generate the data file labeled `lorenz_data.dat`
3. Run the python script `python lorenz_visualization.py`. This visualization will overlay the Fortran data with the Python data and will also show residuals for the state variables.

