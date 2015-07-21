module model
##
 # Copyright (C) 2015 Dynare Team
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

export modeldescription

type modeldescription
    fname::ASCIIString
    dname::ASCIIString
    exonames::Array{ASCIIString,1}
    tex_exonames::Array{ASCIIString,1}
    long_exonames::Array{ASCIIString,1}
    exodetnames::Array{ASCIIString,1}
    tex_exodetnames::Array{ASCIIString,1}
    long_exodetnames::Array{ASCIIString,1}
    endonames::Array{ASCIIString,1}
    tex_endonames::Array{ASCIIString,1}
    long_endonames::Array{ASCIIString,1}
    paramnames::Array{ASCIIString,1}
    tex_paramnames::Array{ASCIIString,1}
    long_paramnames::Array{ASCIIString,1}
    aux_vars::Array{ASCIIString,1}
    exo_nbr::Int
    endo_nbr::Int
    param_nbr::Int
    orig_endo_nbr::Int
    orig_eq_nbr::Int
    eq_nbr::Int
    ramsey_eq_nbr::Int
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
    lead_lag_incidence::Matrix{Int}
    NNZDerivatives::Vector{Int}
    static_and_dynamic_models_differ::Bool
    equations_tags::Array{ASCIIString,1}
    exo_names_orig_ord::Array{Int, 1}
    Sigma_e::Matrix{Float64}
    Correlation_matrix::Matrix{Float64}
    H::Matrix{Float64}
    Correlation_matrix_ME::Matrix{Float64}
    sigma_e_is_diagonal::Bool
    params::Vector{Float64}
end

function modeldescription()
    return modeldescription("",                    # fname
                            "",                    # dname
                            Array(ASCIIString,0),  # exonames
                            Array(ASCIIString,0),  # t_exonames
                            Array(ASCIIString,0),  # l_exonames
                            Array(ASCIIString,0),  # exodetnames
                            Array(ASCIIString,0),  # t_exodetnames
                            Array(ASCIIString,0),  # l_exodetnames
                            Array(ASCIIString,0),  # endonames
                            Array(ASCIIString,0),  # t_endonames
                            Array(ASCIIString,0),  # l_endonames
                            Array(ASCIIString,0),  # paramnames
                            Array(ASCIIString,0),  # t_paramnames
                            Array(ASCIIString,0),  # l_paramnames
                            Array(ASCIIString,0),  # aux_vars
                            0,                     # exo_nbr
                            0,                     # endo_nbr
                            0,                     # param_nbr
                            0,                     # orig_endo_nbr
                            0,                     # orig_eq_nbr
                            0,                     # eq_nbr
                            0,                     # ramsey_eq_nbr
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
                            Array(Int, 3, 0),      # lead_lag_incidence
                            zeros(Int, 3),         # NNZDerivatives
                            false,                 # static_and_dynamic_models_differ
                            Array(ASCIIString,0),  # equations_tags
                            Array(Int64,1),        # exo_names_orig_ord
                            Array(Float64, 0, 0),  # Sigma_e (Covariance matrix of the structural innovations)
                            Array(Float64, 0, 0),  # Correlation_matrix (Correlation matrix of the structural innovations)
                            Array(Float64, 0, 0),  # H (Covariance matrix of the measurement errors)
                            Array(Float64, 0, 0),  # Correlation_matrix_ME (Covariance matrixof the measurement errors)
                            true,                  # sigma_e_is_diagonal
                            Array(Float64, 0)      # params
                            )
end

end
