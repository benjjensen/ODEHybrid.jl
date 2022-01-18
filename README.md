# ODEHybrid.jl

[![CI](https://github.com/benjjensen/ODEHybrid.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benjjensen/ODEHybrid.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/benjjensen/ODEHybrid.jl/branch/main/graph/badge.svg?token=XFGG3UQSNX)](https://codecov.io/gh/benjjensen/ODEHybrid.jl) 
[![](https://img.shields.io/badge/docs-stable-blue.svg)](http://benjjensen.github.io/ODEHybrid.jl/dev)

This is a Julia implementation of the ODE Hybrid package made for MATLAB.
This ported version tries to follow the MATLAB implementation closely so as to maintain the same functionality
and match the MATLAB output 

(NOTE that because the original code does not allow for Events to be used for the custom RK functions, 
relying instead on the ODE suite of solvers in MATLAB like ODE45, Events are currently not supported).

The original package is well-documented at http://www.anuncommonlab.com/doc/odehybrid 