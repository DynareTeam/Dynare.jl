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

The package is tested againt Julia 0.5.x.

## A rational expectations example

The following example is a translation of `example1.mod` that is shipped with
Dynare. It computes the steady state and first order decision rules.

Note that the syntax for the model block has been kept as close as possible to
the original format.

```julia
using Dynare

# Create a model object
m =
@modfile begin
    @var y c k a h b
    @varexo e u
    @parameters beta rho alpha delta theta psi tau
    @model begin
        c*theta*h^(1+psi) = (1-alpha)*y
        k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))*(exp(b(+1))*alpha*y(+1)+(1-delta)*k))
        y = exp(a)*(k(-1)^alpha)*(h^(1-alpha))
        k = exp(b)*(y-c)+(1-delta)*k(-1)
        a = rho*a(-1)+tau*b(-1) + e
        b = tau*a(-1)+rho*b(-1) + u
    end
end

# Do some preliminary computations
compute_model_info(m)

# Define a calibration and some starting values for the nonlinear solver

calib = Dict(
         :alpha => 0.36,
         :rho   => 0.95,
         :tau   => 0.025,
         :beta  => 0.99,
         :delta => 0.025,
         :psi   => 0.0,
         :theta => 2.95
        )

initval = Dict(
           :y => 1.08068253095672,
           :c => 0.80359242014163,
           :h => 0.29175631001732,
           :k => 11.08360443260358,
           :a => 0.0,
           :b => 0.0
           )

# Compute and print the steady state for the given calibration

s = steady_state(m, calib, initval)

print_steady_state(m, s)

# Compute and print eigenvalues and first order decision rules

(gy, gu, eigs) = decision_rules(m, calib, s)

println("Eigenvalues: ", eigs)
println()

print_decision_rules(m, gy, gu)
```

## A perfect foresight example

The following example is a translation of `ramst.mod` that is shipped with
Dynare. It computes the steady state and a perfect foresight simulation.

Note that we had to introduce an endogenous variable `x1` equal to the
exogenous variable `x`, because leads on exogenous variables are not yet supported.

```julia
using Dynare

m =
@modfile begin
    @var c k x1
    @varexo x
    @parameters alph gam delt bet aa
    @model begin
        c + k - aa*x*k(-1)^alph - (1-delt)*k(-1)
        c^(-gam) - (1+bet)^(-1)*(aa*alph*x1(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam)
        x1 = x
    end
end

compute_model_info(m)

calib = Dict(
         :alph => 0.5,
         :gam => 0.5,
         :delt => 0.02,
         :bet => 0.05,
         :aa => 0.5
        )

exoval = Dict( :x => 1.0 )

initval = Dict{Symbol, Float64}()
initval[:k] = ((calib[:delt]+calib[:bet])/(exoval[:x]*calib[:aa]*calib[:alph]))^(1/(calib[:alph]-1))
initval[:c] = calib[:aa]*initval[:k]^calib[:alph]-calib[:delt]*initval[:k]
initval[:x1] = exoval[:x]

s = steady_state(m, calib, initval, exoval)

# Compute a 200 periods perfect foresight simulation
# Start and end at the steady state corresponding to x=1
# There is a shock x=1.2 in t=1
# The results will be in endopath (which includes initial and terminal
#  condition, i.e. it goes from t=0 to t=201)

T = 200
endopath = repmat(s, 1, T+2)
exopath = ones(m.n_exo, T)
exopath[1, 1] = 1.2

perfect_foresight_simul!(m, endopath, exopath, calib)
```
