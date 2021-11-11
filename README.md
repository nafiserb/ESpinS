![alt text](user-guide/images/logo.png)

# ESpinS - Esfahan Spin Simulation package


Package for classical Monte-Carlo simulation of spin systems.
In ESpinS, the spin model hamiltonian of a spin system can be defined through 
the Heisenberg exchange, bi-quadratic, Dzyaloshinskii-Moriya, and single-ion interactions. 
Spin configurations can be updated through either local Metropolis updating or parallel tempering. 

If you are using this software, please cite
> Nafise Rezaei, Mojtaba Alaei, Hadi Akbarzadeh, ESpinS: A program for classical Monte-Carlo simulations of spin systems, [Computational Materials Science, 202, 110947(2022)](https://doi.org/10.1016/j.commatsci.2021.110947)


## Requirements

The package requires the following software:

* Fortran 90 compiler, such as GFortran or IFORT
* BLAS and LAPACK library
* MPI libray (openmpi)

## Installation

    git clone https://github.com/nafiserb/ESpinS.git
    cd ESpinS
    make

In the case that compilation process encounters errors, please edit make.sys file. 
For more information about installation please see user-guide.

## About

Created and maintained by Nafise Rezaei (Physics Department, Isfahan University of Technology, Isfahan, Iran)

With contributions by
* Mojtaba Alaei (Physics Department, Isfahan University of Technology, Isfahan, Iran)
