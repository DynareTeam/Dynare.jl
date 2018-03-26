[![Build Status](https://travis-ci.org/DynareTeam/Dynare.jl.svg?branch=master)](https://travis-ci.org/DynareTeam/Dynare.jl)

## Dynare for Julia

This package aims at bringing to Julia some of the functionality provided by
[Dynare](http://www.dynare.org), a platform for solving economic models and in
particular DSGE models.

Please note that this Julia package is very incomplete compared to the original
Dynare for MATLAB/Octave, but hopefully it will become more featureful over
time.

For the moment the package is only able to compute a model’s steady state,
first order decision rules and perfect foresight simulations.

The package is tested against Julia 0.6.x.

## Installation

Assuming you already have [julia](https://julialang.org/) (version
0.6.x), just do:
```julia-repl
	julia> Pkg.clone("https://github.com/DynareTeam/Dynare.jl", "Dynare")
	julia> Pkg.clone("Dynare", "with-preprocessor")
	julia> Pkg.build("Dynare")
```
The last command will install the preprocessor, selecting the binary
corresponding to your platform (Linux/Windows/MacOS). If you do not
already have them, other packages will be fetched during the installation.

