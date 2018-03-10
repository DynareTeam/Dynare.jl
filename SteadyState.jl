module SteadyState

##
 # Copyright (C) 2016 Dynare Team
 #
 # This file is part of Dynare.
 #
 # Dynare is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # Dynare is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
##

using NLsolve

import DynareModel.Model
import DynareOutput.Output

export steady, steady!
export steady_state, steady_state!

function steady(model::Model, oo::Output)
    if model.analytical_steady_state || model.user_written_analytical_steady_state
        steadystate = zeros(length(model.endo))
        model.steady_state(steadystate, oo.exo_steady_state, model.params)
        return steadystate
    else
        error("You have to provide a closed form solution for the steady state, or declare a guess\nfor the steady state as a third input argument.")
    end
end

function steady!(model::Model, oo::Output)
    if model.analytical_steady_state || model.user_written_analytical_steady_state
        model.steady_state(oo.steady_state, oo.exo_steady_state, model.params)
        return
    else
        error("You have to provide a closed form solution for the steady state, or declare a guess\nfor the steady state as a third input argument.")
    end
end

function steady(model::Model, oo::Output, yinit::Vector{Float64})
    f!(fval::Vector{Float64}, y::Vector{Float64}) = model.static(y, oo.exo_steady_state, model.params, fval)
    j!(fjac::Matrix{Float64}, y::Vector{Float64}) = model.static(y, oo.exo_steady_state, model.params, fjac)
    fj!(fval::Vector{Float64}, fjac::Matrix{Float64}, y::Vector{Float64}) = model.static(y, oo.exo_steady_state, model.params, fval, fjac)
    r = nlsolve(f!, j!, fj!, yinit, show_trace=false)
    if converged(r)
        return r.zero
    else
        return fill(NaN, length(yinit))
    end
end

function steady!(model::Model, oo::Output, yinit::Vector{Float64})
    f!(fval::Vector{Float64}, y::Vector{Float64}) = model.static(y, oo.exo_steady_state, model.params, fval)
    j!(fjac::Matrix{Float64}, y::Vector{Float64}) = model.static(y, oo.exo_steady_state, model.params, fjac)
    fj!(fval::Vector{Float64}, fjac::Matrix{Float64}, y::Vector{Float64}) = model.static(y, oo.exo_steady_state, model.params, fval, fjac)
    r = nlsolve(f!, j!, fj!, yinit, show_trace=false)
    if converged(r)
        oo.steady_state = r.zero
    else
        oo.steady_state = fill(NaN, length(yinit))
    end
end

function steady_state(model::Model, oo::Output)
    ys = steady(model, oo)
    display_steady_state(model, oo, ys)
end

function steady_state!(model::Model, oo::Output)
    steady!(model, oo)
    display_steady_state(model, oo, oo.steady_state)
end

function display_steady_state(model::Model, oo::Output, ys::Vector{Float64})
    println("\n\nSTEADY STATE:\n")
    for i = 1:length(model.endo)
        println(string(model.endo[i].name,  " = ",  ys[i]))
    end
end

function issteadystate(model::Model, oo::Output, ys::Vector{Float64})
    residuals = zeros(Float64, length(ys))
    compute_static_model_residuals!(model, oo, ys, residuals)
    return maximum(abs(residuals))<1e-6
end

function compute_static_model_residuals!(model::Model, oo::Output, ys::Vector{Float64}, residuals::Vector{Float64})
    model.static(ys, oo.exo_steady_state, model.params, residuals)
end

end
