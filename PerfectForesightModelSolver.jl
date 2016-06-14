module PerfectForesightModelSolver

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

import DynareModel.Model
import DynareOutput.Output
import DynareOptions.Options

export simulate_perfect_foresight_model!

function simulate_perfect_foresight_model!(endogenousvariables::Matrix{Float64}, exogenousvariables::Matrix{Float64}, steadystate::Vector{Float64}, model::Model, options::Options)

    lead_lag_incidence = model.lead_lag_incidence

    nyp = countnz(lead_lag_incidence[1,:])
    ny0 = countnz(lead_lag_incidence[2,:])
    nyf = countnz(lead_lag_incidence[3,:])

    ny = length(model.endo)
    nd = nyp+ny0+nyf

    periods = options.pfmsolver.periods
    params = model.params

    tmp = lead_lag_incidence[2:3,:]'
    i_cols_A1 = find(tmp)
    i_cols_1  = tmp[i_cols_A1]

    tmp = lead_lag_incidence[1:2,:]'
    i_cols_AT = find(tmp)
    i_cols_T  = tmp[i_cols_AT]

    tmp = lead_lag_incidence[2,:]'
    i_cols_A0 = find(tmp)
    i_cols_0  = tmp[i_cols_A0]

    i_cols_j = collect(1:nd)
    i_upd = ny+collect(1:periods*ny)

    Y = vec(endogenousvariables)
    z = Y[find(lead_lag_incidence')]

    jacobian = zeros(Float64, ny, length(z)+length(model.exo))
    residuals = zeros(Float64, ny)

    println("\nMODEL SIMULATION:\n")

    rd = zeros(Float64, periods*ny)

    iA = zeros(Int64, periods*model.nnzderivatives[1])
    jA = zeros(Int64, periods*model.nnzderivatives[1])
    vA = zeros(Float64, periods*model.nnzderivatives[1])

    convergence = false
    iteration = 0
    
    while !convergence
        iteration += 1
        i_rows = collect(1:ny)
        i_cols_A = find(lead_lag_incidence')
        i_cols = i_cols_A
        m = 0
        for it = 2:(periods+1)
            model.dynamic(Y[i_cols], exogenousvariables, params, steadystate, it, residuals, jacobian)
            if it==(periods+1) & it==2
                (r, c, v) = findnz(jacobian[:,i_cols_0])
                k = collect(1:length(v))+m
                iA[k] = i_rows[r]
                jA[k] = i_cols_A0[c]
                vA[k] = v
            elseif it==(periods+1)
                (r, c, v) = findnz(jacobian[:,i_cols_T])
                k = collect(1:length(v))+m
                iA[k] = i_rows[r]
                jA[k] = i_cols_A[i_cols_T[c]]
                vA[k] = v
            elseif it==2
                (r, c, v) = findnz(jacobian[:,i_cols_1])
                k = collect(1:length(v))+m
                iA[k] = i_rows[r]
                jA[k] = i_cols_A1[c]
                vA[k] = v
            else
                (r, c, v) = findnz(jacobian[:,i_cols_j])
                k = collect(1:length(v))+m
                iA[k] = i_rows[r]
                jA[k] = i_cols_A[c]
                vA[k] = v
            end
            m += length(v)
            rd[i_rows] = residuals
            i_rows += ny
            i_cols += ny
            if it>2
                i_cols_A += ny
            end
        end
        err = maximum(abs(rd))
        println("Iter. ", iteration, "\t err. ", round(err, 12))
        if err<options.pfmsolver.tolf
            iteration -= 1
            convergence = true
        end
        A = sparse(iA[1:m], jA[1:m], vA[1:m])
        dy = -A\rd
        Y[i_upd] += dy
        if maximum(abs(dy))<options.pfmsolver.tolx
            convergence = true
        end
    end
    if convergence
        println("\nPFM solver converged in ", iteration, " iterations!\n")
        endogenousvariables = reshape(Y, ny, periods+2)
    end
end

end
