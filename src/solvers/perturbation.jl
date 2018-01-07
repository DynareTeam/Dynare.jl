using ..DynLinAlg.LinSolveAlgo
using model

struct ResultsPerturbationWs
    g::Array{Matrix{Float64}}  # full approximation
    gs::Array{Matrix{Float64}} # state transition matrices
    g1_1::SubArray # solution first order derivatives w.r. to state variables
    g1_2::SubArray # solution first order derivatives w.r. to current exogenous variables
    g1_3::SubArray # solution first order derivatives w.r. to lagged exogenous variables
    f1g1plusf2::Matrix{Float64} # f_1*g_1 + f_2
    f1g1plusf2_linsolve_ws::LinSolveWS # LU decomposition of f_1*g_2 + f_2
    
    function ResultsPerturbationWs(m::Model,order::Int64)
        nstate = m.n_bkwrd + m.n_both
        g =  [Matrix{Float64}(m.endo_nbr,(m.n_bkwrd + m.n_both + m.current_exogenous_nbr + 1)^k) for k = 1:order]
        gs = [Matrix{Float64}(m.n_bkwrd+m.n_both,(m.n_bkwrd+m.n_both)^k) for k = 1:order]
        g1_1 = view(g[1],:,1:nstate)
        g1_2 = view(g[1],:,nstate + (1:m.current_exogenous_nbr))
        g1_3 = view(g[1],:,nstate + m.current_exogenous_nbr + (1:m.lagged_exogenous_nbr))
        f1g1plusf2 = Matrix{Float64}(m.endo_nbr,m.endo_nbr)
        f1g1plusf2_linsolve_ws = LinSolveWS(m.endo_nbr)
        new(g, gs, g1_1, g1_2, g1_3, f1g1plusf2, f1g1plusf2_linsolve_ws)
    end
end

