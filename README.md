[![DOI](https://zenodo.org/badge/375122337.svg)](https://zenodo.org/badge/latestdoi/375122337)

# Phase reduction analysis of periodic thermoacoustic oscillations in a Rijke tube (companion code)

This repository contains companion code for the article ["Phase reduction analysis of periodic thermoacoustic oscillations in a Rijke tube"](https://www.doi.org/10.1017/jfm.2021.1093) by C. S. Skene, K. Taira (JFM, 2022).

The code is written in Julia and utilises the following libraries

* Conda
* IJulia
* SparseArrays
* LinearAlgebra
* QuadGK
* DelayDiffEq
* NonlinearEigenproblems
* Statistics
* FileIO
* JLD
* JLD2
* PyPlot
* FFTW

The file _GalerkinFuncs.jl_, contains the structure needed to create the Galerkin model for a given set of parameters as well as useful functions for the phase reduction analysis. The initial conditions used for the paper are contained in the folder _data_. Also provided are Jupyter notebooks that use this file to perform the following tasks.

## 1-Neutral Curve
The notebook _1-NeutralCurve.ipynb_ setups the linearised equations and solves an eigenvalue problem in order to find the neutral curve.

## 2-Phase Sensitivity
The notebook _2-PhaseSensitivity.ipynb_ solves the non-linear equations to find the limit cycle. Linearising about this periodic solution, the adjoint equations are then solved in order to find the phase sensitivity function. The correct normalisation is then found via the bilinear form.

## 3-Phase Coupling Function
The notebook _3-PhaseCouplingFunction.ipynb_ loads the saved phase sensitivity solution from notebook _2-PhaseSensitivity.ipynb_ and finds the phase coupling function for the global forcing considered in the paper. Using this phase coupling function the Arnold tongues for m:n phase locking are found.

## 4-Parametric Sensitivity
The notebook _4-ParametricSensitivity.ipynb_ computes the synchronisability as the location of the pressure actuation in the Rijke tube is moved. The procedure used in notebooks 2 and 3 is automated in order to carry out this procedure for a range of parameters.

# Citation
If you find this repository useful for your research please cite the paper
```
@article{skene_taira_2022,
title={Phase-reduction analysis of periodic thermoacoustic oscillations in a Rijke tube},
volume={933},
DOI={10.1017/jfm.2021.1093},
journal={Journal of Fluid Mechanics},
publisher={Cambridge University Press},
author={Skene, Calum S. and Taira, Kunihiko},
year={2022},
pages={A35}}
```
# Acknowledgements
This work was supported by the US Air Force Office of Scientific Research (FA9550-16-1-0650 and FA9550-21-1-0178, monitored by Douglas Smith and Gregg Abate).
