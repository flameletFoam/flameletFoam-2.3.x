flameletFoam-2.3.x
==================

## Description

This code realizes a steady laminar flamelet approach for turbulent non-premixed combustion.
The solver is based on ''rhoReactingFoam'', i.e. it is pressure-based (PISO), compressible and runs with both LES and RAS turbulence.
 
The theory is mainly taken from the work of N. Peters and is based on the view of a turbulent flame as an ensemble of laminar flamelets.
The calculation of these flamelets is a one-dimensional problem and can be done in a pre-processing step.
Integration using a presumed beta-Probability Density Function (PDF) accounts for the interaction between turbulent fluctuations and chemistry.
The results of the pre-processing procedure are stored in tables which are accessed during the simulation.
Values of interest, such as species mass fraction or enthalpy, are looked-up and connected to the flow using three parameters - the mixture fraction, its variance and the scalar dissipation rate.
In doing so, the expensive solution of chemical mechanisms during run-time can be avoided and the run-time thus reduces significantly.

More information is available on the Extend-bazaar page:
https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/flameletFoam

## Installation

This version works with OpenFOAM-2.3

* Prepare a directory on your system, e.g.:  

  `mkdir ~/OpenFOAM/flamletFoam/`

* Download flameletFoam using git:

  `git clone https://github.com/flameletFoam/flameletFoam-2.3.x/ ~/OpenFOAM/flameletFoam/`

* Set an environment variable to the flameletFoam src folder:

  `export LIB_FLAMELET_SRC=~/OpenFOAM/flameletFoam/src/`

* Execute `./Allwmake`

## Tutorials

There is a RANS and a LES tutorial available:

  `cd ./tutorials/pilotedDiffusionFlame/ras/`

  `./Allrun`

## Notes

More information on using flameletFoam, in particular a description of the workflow including the table-generation with cantera, is available on:
https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/flameletFoam

Please feel free to contact me should you find bugs or have suggestions how to make the code better.
