module FirstOrderSolver

import ...model: Model, get_de, get_abc
using ...DynLinAlg.QrAlgo
using ..Solvers.CyclicReduction
using ..Solvers.GeneralizedSchurDecompositionSolver
import ...DynLinAlg.LinSolveAlgo: LinSolveWS, linsolve_core!, linsolve_core_no_lu!, lu!
import ..Solvers: ResultsPerturbationWs
import ..SolveEyePlusMinusAkronB: EyePlusAtKronBWS, generalized_sylvester_solver!

import Base.LinAlg.BLAS: gemm!

export FirstOrderSolverWS, first_order_solver

type FirstOrderSolverWS
    jacobian_static::Matrix{Float64} 
    qr_ws::QrWS
    solver_ws::Union{GsSolverWS,CyclicReductionWS}
    ghx::StridedMatrix{Float64}
    gx::Matrix{Float64}
    hx::Matrix{Float64}
    temp1::Matrix{Float64}
    temp2::Matrix{Float64}
    temp3::Matrix{Float64}
    temp4::Matrix{Float64}
    temp5::Matrix{Float64}
    b10::Matrix{Float64}
    b11::Matrix{Float64}
    linsolve_static_ws::LinSolveWS
    eye_plus_at_kron_b_ws::EyePlusAtKronBWS
    
    function FirstOrderSolverWS(algo::String, jacobian::Matrix, m::Model)
        if m.n_static > 0
            jacobian_static = Matrix{Float64}(m.endo_nbr,m.n_static)
            qr_ws = QrWS(jacobian_static)
        else
            jacobian_static = Matrix{Float64}(0,0)
            qr_ws = QrWS(Matrix{Float64}(0,0))
        end
        if algo == "GS"
            d = zeros(m.n_dyn,m.n_dyn)
            e = zeros(m.n_dyn,m.n_dyn)
            solver_ws = GsSolverWS(d,e,m.n_bkwrd+m.n_both)
        elseif algo == "CR"
            n = m.endo_nbr - m.n_static
            solver_ws = CyclicReductionWS(n)
        end
        ghx = Matrix{Float64}(m.endo_nbr,m.n_bkwrd+m.n_both)
        gx = Matrix{Float64}(m.n_fwrd+m.n_both,m.n_bkwrd+m.n_both)
        hx = Matrix{Float64}(m.n_bkwrd+m.n_both,m.n_bkwrd+m.n_both)
        temp1 = Matrix{Float64}(m.n_static,m.n_fwrd+m.n_both)
        temp2 = Matrix{Float64}(m.n_static,m.n_bkwrd+m.n_both)
        temp3 = Matrix{Float64}(m.n_static,m.n_bkwrd+m.n_both)
        temp4 = Matrix{Float64}(m.endo_nbr - m.n_static,m.n_bkwrd+m.n_both)
        temp5 = Matrix{Float64}(m.endo_nbr,max(m.current_exogenous_nbr,m.lagged_exogenous_nbr))
        b10 = Matrix{Float64}(m.n_static,m.n_static)
        b11 = Matrix{Float64}(m.n_static,length(m.p_current_ns))
        linsolve_static_ws = LinSolveWS(m.n_static)
        if m.serially_correlated_exogenous
            eye_plus_at_kron_b_ws = EyePlusAtKronBWS(ma, mb, mc, 1)
        else
            eye_plus_at_kron_b_ws = EyePlusAtKronBWS(1, 1, 1, 1)
        end
        new(jacobian_static, qr_ws, solver_ws, ghx, gx, hx, temp1,
            temp2, temp3, temp4, temp5, b10, b11, linsolve_static_ws,
            eye_plus_at_kron_b_ws)
    end
end
        
function remove_static!(ws::FirstOrderSolverWS,jacobian::Matrix,p_static::Vector)
    ws.jacobian_static[:,:] = view(jacobian,:,p_static)
    dgeqrf_core!(ws.qr_ws,ws.jacobian_static)
    dormrqf_core!(ws.qr_ws,Ref{UInt8}('L'),Ref{UInt8}('T'),ws.jacobian_static,
                  jacobian)
end

function add_static!(results::ResultsPerturbationWs,ws::FirstOrderSolverWS,jacobian::Matrix{Float64},model::Model)
    i_static = 1:model.n_static
    #    temp = - jacobian[i_static,model.p_fwrd_b]*gx*hx
    ws.temp1 = view(jacobian,i_static,model.p_fwrd_b)
    gemm!('N','N',1.0,ws.temp1,ws.gx,0.0,ws.temp2)
    gemm!('N','N',-1.0,ws.temp2,ws.hx,0.0,ws.temp3)
    ws.b10 = view(jacobian,i_static, model.p_static)
    ws.b11 = view(jacobian,i_static, model.p_current_ns)
    ws.temp3 -= view(jacobian,i_static,model.p_bkwrd_b)
    for i=1:(model.n_bkwrd + model.n_both)
        for j=1:length(model.i_dyn)
            ws.temp4[j,i] = results.g[1][model.i_dyn[j],i]
        end
    end
    gemm!('N','N',-1.0,ws.b11,ws.temp4,1.0,ws.temp3)
    linsolve_core!(ws.linsolve_static_ws,Ref{UInt8}('N'),ws.b10,ws.temp3)
    for i = 1:model.n_states
        for j=1:model.n_static
            results.g[1][model.i_static[j],i] = ws.temp3[j,i]
        end
    end
end

using Base.Test
function make_f1g1plusf2!(results::ResultsPerturbationWs,model,jacobian)
    nstate = model.n_bkwrd + model.n_both
    so = nstate*model.endo_nbr + 1
    for i=1:model.n_current
        copy!(results.f1g1plusf2,(model.i_current[i]-1)*model.endo_nbr+1,jacobian,so,model.endo_nbr)
        so += model.endo_nbr
    end
    offset = nstate + model.n_current
    for i = 1:nstate
        y = view(results.g[1],:,i)
        z = view(results.f1g1plusf2, :, model.i_bkwrd_b[i])
        for j=1:nstate
            x = 0.0
            for k=1:(model.n_fwrd + model.n_both)
                x += jacobian[j, offset + k]*y[model.i_fwrd_b[k]]
            end
            z[j] += x
        end
    end
    lu!(results.f1g1plusf2_linsolve_ws.lu, results.f1g1plusf2, results.f1g1plusf2_linsolve_ws.ipiv)
end

function solve_for_derivatives_with_respect_to_shocks(results::ResultsPerturbationWs, jacobian::AbstractMatrix, ws::FirstOrderSolverWS, model::Model)
    if model.lagged_exogenous_nbr > 0
        f6 = view(jacobian,:,model.i_lagged_exogenous)
        for i = 1:model.current_exogenous_nbr
            for j = 1:model.endo_nbr
                results.g1_3[i,j] = -f6[i,j]
            end
        end
        linsolve_core_no_lu!(ws, Ref{UInt8}('N'), results.f1g1plusf2, results.g1_3)
    end
    if model.current_exogenous_nbr > 0
        copy!(results.g[1], (model.n_bkwrd + model.n_both)*model.endo_nbr + 1, jacobian,
              (model.n_fwrd + model.n_bkwrd + 2*model.n_both + model.n_current)*model.endo_nbr + 1,
              model.current_exogenous_nbr*model.endo_nbr)
        gu = view(results.g[1],:, model.n_bkwrd + model.n_both + (1:model.current_exogenous_nbr))
        scale!(gu, -1.0)
        if model.serially_correlated_exogenous
            # TO BE DONE
        else
            linsolve_core_no_lu!(results.f1g1plusf2_linsolve_ws, Ref{UInt8}('N'), results.f1g1plusf2, gu)
        end
    end
end

function first_order_solver(results::ResultsPerturbationWs,ws::FirstOrderSolverWS,algo::String, jacobian::Matrix, model::Model, options)
    if model.n_static > 0
        remove_static!(ws,jacobian,model.p_static)
    end
    n = model.n_fwrd + model.n_bkwrd + model.n_both
    if algo == "CR"
        a, b, c = get_abc(model,jacobian)
	x = zeros(a)
        cyclic_reduction!(x,a,b,c,ws.solver_ws,options.cycle_reduction.tol,100)
        if ws.solver_ws.info[1] > 0
            error("CR didn't converge")
        end
        for i = 1:length(model.i_bkwrd_ns)
            for j = 1:length(model.i_dyn)
                results.g[1][model.i_dyn[j],i] = x[j,model.i_bkwrd_ns[i]]
            end
            for j = 1:length(model.i_bkwrd_ns)
                results.gs[1][j,i] = results.g[1][model.hx_rows[j],i]
            end
        end
        
    elseif algo == "GS"
        d, e = get_de(jacobian[model.n_static+1:end,:],model)
        gs_solver!(ws.solver_ws,d,e,model.n_bkwrd+model.n_both,options.generalized_schur.criterium)
        results.gs[1] = ws.solver_ws.g2
        for i = 1:model.n_bkwrd+model.n_both
            for j = 1:model.n_bkwrd
                results.g[1][model.i_bkwrd[j],i] = ws.solver_ws.g1[j,i]
            end
            for j = 1:model.n_fwrd
                results.g[1][model.i_fwrd[j],i] = ws.solver_ws.g2[j,i]
            end
            for j = 1:model.n_both
                results.g[1][model.i_both[j],i] = ws.solver_ws.g2[model.n_fwrd+j,i]
            end
        end
    end
    if model.n_static > 0
        add_static!(results, ws, jacobian, model)
    end
    make_f1g1plusf2!(results, model, jacobian)        
    solve_for_derivatives_with_respect_to_shocks(results, jacobian, ws, model)
end

end    
