module SteadyState

##
 # Copyright (C) 2016-2018 Dynare Team
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

using DynareSolvers
using AutoAligns

import DynareModel.Model
import DynareOutput.Output

export staticmodule
export steady, steady!
export steady_state, steady_state!

function steady(model::Model, oo::Output)
    if model.analytical_steady_state || model.user_written_analytical_steady_state
        steadystate = zeros(model.endo_nbr)
        model.steady_state(steadystate, oo.exo_steady_state, model.params)
        return steadystate
    else
        error("You have to provide a closed form solution for the steady state, or declare a guess\nfor the steady state as a third input argument.")
    end
end

function steady!(model::Model, oo::Output)
    if model.analytical_steady_state || model.user_written_analytical_steady_state
        model.steady_state(oo.steady_state, oo.exo_steady_state, model.params)
    else
        error("You have to provide a closed form solution for the steady state, or declare a guess\nfor the steady state as a third input argument.")
    end
end

function steady(model::Model, oo::Output, yinit::Vector{Float64})
    T = zeros(Float64, sum(model.temporaries.static[1:2]))
    f!(fval::Vector{Float64}, y::Vector{Float64}) = model.static(T, fval, y , oo.exo_steady_state, model.params)
    j!(fjac::Matrix{Float64}, y::Vector{Float64}) = model.static(T, fjac, y, oo.exo_steady_state, model.params)
    ys, info = trustregion(f!, j!, yinit, 1.0e-6, 1.0e-6, 100)
    return ys, info
end

function steady!(model::Model, oo::Output, yinit::Vector{Float64})
    T = zeros(Float64, sum(model.temporaries.static[1:2]))
    f!(fval::Vector{Float64}, y::Vector{Float64}) = model.static(T, fval, y , oo.exo_steady_state, model.params)
    j!(fjac::Matrix{Float64}, y::Vector{Float64}) = model.static(T, fjac, y, oo.exo_steady_state, model.params)
    oo.steady_state, info = trustregion(f!, j!, yinit, 1.0e-6, 1.0e-6, 100)
end

function steady_state(model::Model, oo::Output, yinit::Vector{Float64}, display::Bool=true)
    ys, info = steady(model, oo, yinit)
    if info==1
        if display
            display_steady_state(model, oo, ys)
        end
    else
        error("Steady state not found!")
    end
    return ys, info
end

function steady_state!(model::Model, oo::Output, yinit::Vector{Float64}, display::Bool=true)
    steady!(model, oo, yinit)
    if issteadystate(model, oo, oo.steady_state)
        if display
            display_steady_state(model, oo)
        end
    else
        error("Steady state not found!")
    end
end

function steady_state(model::Model, oo::Output, display::Bool=true)
    ys = steady(model, oo)
    if !issteadystate(model, oo, ys)
        error("Steady state provided in steady state block/file is not correct!")
    end
    if display
        display_steady_state(model, oo, ys)
    end
    return ys, info
end

function steady_state!(model::Model, oo::Output, display::Bool=true)
    steady!(model, oo)
    if !issteadystate(model, oo, oo.steady_state)
        error("Steady state provided in steady state block/file is not correct!")
    end
    if issteadystate(oo.steady_state)
        if display
            display_steady_state(model, oo)
        end
    else
        error("Steady state not found!")
    end
end

function display_steady_state(model::Model, oo::Output)
    aa = AutoAlign(align = Dict(1 => left, 2=> left, :default => right))
    println("\n  Deterministic Steady State:\n")
    for i = 1:model.endo_nbr
        if abs(oo.steady_state[i])<1.0e-6
            tmp = zero(Float64)
        else
            tmp = oo.steady_state[i]
        end
        print(aa, "    ", model.endo[i].name, "\t=\t", @sprintf("%.6f", tmp))
        println(aa)
    end
    print(STDOUT, aa)
end

function display_steady_state(model::Model, oo::Output, ys::Vector{Float64})
    aa = AutoAlign(align = Dict(1 => left, 2=> left, :default => right))
    println("\n  Deterministic Steady State:\n")
    for i = 1:model.endo_nbr
        if abs(ys[i])<1.0e-6
            tmp = zero(Float64)
        else
            tmp = ys[i]
        end
        print(aa, "    ", model.endo[i].name, "\t=\t", @sprintf("%.6f", tmp))
        println(aa)
    end
    print(STDOUT, aa)
end

function issteadystate(model::Model, oo::Output, ys::Vector{Float64})
    r = zeros(Float64, model.endo_nbr)
    t = zeros(Float64, model.temporaries.static[1])
    model.static(t, r, ys, oo.exo_steady_state, model.params)
    return norm(r)<1e-6
end

end
