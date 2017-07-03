module FirstOrder

import ...model: Model, get_de, get_abc
using ..LinAlg.qr_algo
using ..CyclicReduction
using ..gs_solver
import ..LinAlg.linsolve_algo: LinSolveWS, linsolve_core!

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
    b10::Matrix{Float64}
    b11::Matrix{Float64}
    linsolve_ws::LinSolveWS
    
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
        b10 = Matrix{Float64}(m.n_static,m.n_static)
        b11 = Matrix{Float64}(m.n_static,length(m.p_current_ns))
        linsolve_ws = LinSolveWS(m.n_static)
        
        new(jacobian_static,qr_ws,solver_ws,ghx,gx,hx,temp1,temp2,temp3,temp4,b10,b11,linsolve_ws)
    end
end
        
function remove_static!(ws::FirstOrderSolverWS,jacobian::Matrix,p_static::Vector)
    ws.jacobian_static[:,:] = view(jacobian,:,p_static)
    dgeqrf_core!(ws.qr_ws,ws.jacobian_static)
    dormrqf_core!(ws.qr_ws,Ref{UInt8}('L'),Ref{UInt8}('T'),ws.jacobian_static,
                  jacobian)
end

function add_static!(ws::FirstOrderSolverWS,jacobian::Matrix{Float64},model::Model)
    i_static = 1:model.n_static
    #    temp = - jacobian[i_static,model.p_fwrd_b]*gx*hx
    ws.temp1 = view(jacobian,i_static,model.p_fwrd_b)
    gemm!('N','N',1.0,ws.temp1,ws.gx,0.0,ws.temp2)
    gemm!('N','N',-1.0,ws.temp2,ws.hx,0.0,ws.temp3)
    ws.b10 = view(jacobian,i_static, model.p_static)
    ws.b11 = view(jacobian,i_static, model.p_current_ns)
    ws.temp3 -= view(jacobian,i_static,model.p_bkwrd_b)
    for i=1:size(ws.ghx,2)
        for j=1:length(model.i_dyn)
            ws.temp4[j,i] = ws.ghx[model.i_dyn[j],i]
        end
    end
    gemm!('N','N',-1.0,ws.b11,ws.temp4,1.0,ws.temp3)
    linsolve_core!(ws.linsolve_ws,Ref{UInt8}('N'),ws.b10,ws.temp3)
    ws.ghx[model.i_static,:] = ws.temp3
end

function first_order_solver(ws::FirstOrderSolverWS,algo::String, jacobian::Matrix, model::Model, options)
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
        ws.ghx[model.i_dyn,:] = view(x,:,model.i_bkwrd_ns)
        ws.gx[:,:] = view(ws.ghx,model.gx_rows,:)
        ws.hx[:,:] = view(ws.ghx,model.hx_rows,:)
    elseif algo == "GS"
        d, e = get_de(jacobian[model.n_static+1:end,:],model)
        gs_solver!(ws.solver_ws,d,e,model.n_bkwrd+model.n_both,options.generalized_schur.criterium)
        ws.gx = ws.solver_ws.g2
        ws.hx = ws.solver_ws.g1
        ws.ghx[model.i_bkwrd,:] = ws.solver_ws.g1[1:model.n_bkwrd,:]
        ws.ghx[model.i_fwrd,:] = ws.solver_ws.g2[1:model.n_fwrd,:]
        ws.ghx[model.i_both,:] = ws.solver_ws.g2[model.n_fwrd+(1:model.n_both),:]
    end
    if model.n_static > 0
        add_static!(ws,jacobian,model)
    end
end

end    
