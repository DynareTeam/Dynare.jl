module DynareModel
##
 # Copyright (C) 2015-2018 Dynare Team
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

export Model, Endo, Exo, ExoDet, Param, dynare_model

abstract type Atom end

immutable Endo <: Atom
    name::String
    tex_name::String
    long_name::String
end

immutable Exo <: Atom
    name::String
    tex_name::String
    long_name::String
end

immutable ExoDet <: Atom
    name::String
    tex_name::String
    long_name::String
end

immutable Param <: Atom
    name::String
    tex_name::String
    long_name::String
end

immutable AuxVars
    endo_index::Int
    var_type::Int
    orig_index::Int
    orig_lead_lag::Int
    eq_nbr::Int
    orig_expr::String
end

immutable PredVars
    index::Int
end

immutable ObsVars
    index::Int
end

immutable DetShocks
    exo_det::Int
    exo_id::Int
    multiplicative::Bool
    periods::Vector{Int}
    value::Float64
end

immutable EquationTag
    eq_nbr::Int
    name::String
    value::String
end

type Model
    fname::String
    dname::String
    dynare_version::String
    endo::Vector{Endo}
    exo::Vector{Exo}
    exo_det::Vector{ExoDet}
    param::Vector{Param}
    aux_vars::Vector{AuxVars}
    pred_vars::Vector{Int}
    obs_vars::Vector{Int}
    state_var::Vector{Int}
    orig_endo_nbr::Int
    orig_eq_nbr::Int
    eq_nbr::Int
    ramsey_eq_nbr::Int
    det_shocks::Vector{DetShocks}
    nstatic::Int
    nfwrd::Int
    npred::Int
    nboth::Int
    nsfwrd::Int
    nspred::Int
    ndynamic::Int
    maximum_lag::Int
    maximum_lead::Int
    maximum_endo_lag::Int
    maximum_endo_lead::Int
    maximum_exo_lag::Int
    maximum_exo_lead::Int
    orig_maximum_lag::Int
    orig_maximum_lead::Int
    orig_maximum_endo_lag::Int
    orig_maximum_endo_lead::Int
    orig_maximum_exo_lag::Int
    orig_maximum_exo_lead::Int
    orig_maximum_exo_det_lag::Int
    orig_maximum_exo_det_lead::Int
    lead_lag_incidence::Matrix{Int}
    nnzderivatives::Vector{Int}
    analytical_steady_state::Bool
    user_written_analytical_steady_state::Bool
    static_and_dynamic_models_differ::Bool
    equation_tags::Vector{String}
    exo_names_orig_ord::Vector{Int}
    sigma_e::Matrix{Float64}
    correlation_matrix::Matrix{Float64}
    h::Matrix{Float64}
    correlation_matrix_me::Matrix{Float64}
    sigma_e_is_diagonal::Bool
    params::Vector{Float64}
    static::Function
    static_params_derivs::Function
    dynamic::Function
    dynamic_params_derivs::Function
    steady_state::Function
end

function dynare_model()
    return Model("",                    # fname
                 "",                    # dname
                 "",                    # dynare_version
                 Vector{Endo}(),        # endo
                 Vector{Exo}(),         # exo
                 Vector{ExoDet}(),      # exo_det
                 Vector{Param}(),       # param
                 Vector{AuxVars}(),     # aux_vars
                 Vector{Int}(),         # pred_vars
                 Vector{Int}(),         # obs_vars
                 Vector{Int}(),         # state_var
                 0,                     # orig_endo_nbr
                 0,                     # orig_eq_nbr
                 0,                     # eq_nbr
                 0,                     # ramsey_eq_nbr
                 Vector{DetShocks}(),   # det_shocks
                 0,                     # nstatic
                 0,                     # nfwrd
                 0,                     # npred
                 0,                     # nboth
                 0,                     # nsfwrd
                 0,                     # nspred
                 0,                     # ndynamic
                 0,                     # maximum_lag
                 0,                     # maximum_lead
                 0,                     # maximum_endo_lag
                 0,                     # maximum_endo_lead
                 0,                     # maximum_exo_lag
                 0,                     # maximum_exo_lead
                 0,                     # orig_maximum_lag
                 0,                     # orig_maximum_lead
                 0,                     # orig_maximum_endo_lag
                 0,                     # orig_maximum_endo_lead
                 0,                     # orig_maximum_exo_lag
                 0,                     # orig_maximum_exo_lead
                 0,                     # orig_maximum_exo_det_lag
                 0,                     # orig_maximum_exo_det_lead
                 Matrix{Int}(0,0),         # lead_lag_incidence
                 zeros(Int64,3),         # nnzderivatives
                 false,                 # analytical_steady_state
                 false,                 # user_written_analytical_steady_state
                 false,                 # static_and_dynamic_models_differ
                 Vector{String}(),      # equation_tags
                 Vector{Int}(),         # exo_names_orig_ord
                 Matrix{Float64}(0,0),     # sigma_e (Cov matrix of the structural innovations)
                 Matrix{Float64}(0,0),     # correlation_matrix (Corr matrix of the structural innovations)
                 Matrix{Float64}(0,0),     # h (Cov matrix of the measurement errors)
                 Matrix{Float64}(0,0),     # correlation_matrix_me (Cov matrix of the measurement errors)
                 true,                  # sigma_e_is_diagonal
                 Vector{Float64}(),     # params
                 function()end,         # static
                 function()end,         # static_params_derivs
                 function()end,         # dynamic
                 function()end,         # dynamic_params_derivs
                 function()end          # steady_state
                )
end

end
