module SecondOrder

import Base.LinAlg.BLAS: gemm!        
import ...model: Model
using ..FirstOrder
using ..GeneralizedSylvester
import ..LinAlg.Kronecker: a_mul_kron_b!

export SecondOrderSolverWS, second_order_solver!

type SecondOrderSolverWS
    a::Matrix{Float64}
    atmp::Matrix{Float64}
    b::Matrix{Float64}
    c::Matrix{Float64}
    d::Matrix{Float64}
    ff::Matrix{Float64}
    fplus::Matrix{Float64}
    ghx_fwrd::Matrix{Float64}
    zz::Matrix{Float64}
    gxx::Matrix{Float64}
    gxu::Matrix{Float64}
    guu::Matrix{Float64}
    gss::Vector{Float64}
    kstates::SubArray{Int64}
    work::Vector{Float64}
    model_derivatives_1_current::SubArray{Float64}
    model_derivatives_1_bkwrd::SubArray{Float64}
    model_derivatives_1_fwrd::SubArray{Float64}
    model_derivatives_1_exo::SubArray{Float64}
    model_derivatives_2_endo::SubArray{Float64}
    zz_fwrd::SubArray{Float64}
    zz_current::SubArray{Float64}
    zz_bkwrd::SubArray{Float64}
    function SecondOrderSolverWS(endo_nbr::Int64, exo_nbr::Int64, model::Model)
        n1 = model.n_bkwrd + model.endo_nbr + model.n_fwrd + 2*model.n_both
        n2 = n1 + exo_nbr
        nstates = model.n_bkwrd + model.n_both
        a = Matrix{Float64}(endo_nbr,endo_nbr)
        atmp = Matrix{Float64}(endo_nbr,nstates)
        b = Matrix{Float64}(endo_nbr,endo_nbr)
        c = Matrix{Float64}(nstates,nstates)
        d = Matrix{Float64}(endo_nbr,nstates*nstates)
        fplus = Matrix{Float64}(endo_nbr,model.n_fwrd + model.n_both)
        ghx_fwrd = Matrix{Float64}(model.n_fwrd + model.n_both, model.n_bkwrd + model.n_both)
        ff = Matrix{Float64}(endo_nbr,n1*n1)
        zz = Matrix{Float64}(n1,nstates)
        gxx = Matrix{Float64}(endo_nbr,nstates*nstates)
        gxu = Matrix{Float64}(endo_nbr,nstates*exo_nbr)
        guu = Matrix{Float64}(endo_nbr,exo_nbr*exo_nbr)
        gss = Vector{Float64}(endo_nbr)
        kk = reshape(1:n2*n2,n2,n2)
        kstates = view(kk,1:n1,1:n1)
        work = Vector{Float64}(endo_nbr*n1*nstates)
        new(a, atmp, b, c, d, ff, fplus, ghx_fwrd, zz, gxx, gxu, guu, gss, kstates, work)
    end
end

function second_order_solver_setup!(model_derivatives_1::AbstractMatrix, model_derivatives_2::AbstractMatrix, model::Model, ws::SecondOrderSolverWS)
    n1 = model.n_bkwrd + model.endo_nbr + model.n_fwrd + 2*model.n_both
#    n2 = n1 + exo_nbr
    ws.model_derivatives_1_current = view(model_derivatives_1, :, model.p_current)
    ws.model_derivatives_1_bkwrd = view(model_derivatives_1, :, model.p_bkwrd_b)
    ws.model_derivatives_1_fwrd = view(model_derivatives_1, :, model.p_fwrd_b)
#    ws.model_derivatives_1_exo = view(model_derivatives_1, :, n1+1:n2)
    ws.model_derivatives_2_endo = view(model_derivatives_2,:,vec(ws.kstates))
    ws.zz_bkwrd =  view(ws.zz,1:model.n_bkwrd + model.n_both,:)
    ws.zz_current = view(ws.zz,model.n_bkwrd+model.n_both+(1:model.n_current),:)
    ws.zz_fwrd = view(ws.zz,model.n_bkwrd+model.n_both+model.n_current+(1:model.n_fwrd+model.n_both),:)
end

function second_order_solver!(model_derivatives_1::AbstractMatrix, model_derivatives_2::AbstractMatrix, model::Model, first_order_ws::FirstOrderSolverWS, ws::SecondOrderSolverWS)
    second_order_solver_setup!(model_derivatives_1, model_derivatives_2, model, ws)
    fill!(ws.a,0.0)
    copy!(ws.a,ws.model_derivatives_1_current)

    copy!(ws.ghx_fwrd,view(first_order_ws.ghx, model.i_fwrd_b, :))
    copy!(ws.fplus,ws.model_derivatives_1_fwrd)
    A_mul_B!(ws.atmp,ws.fplus, ws.ghx_fwrd)
    println("ws.a")
    display(ws.a)
    println("ws.atmp")
    display(ws.atmp)
    for j = 1:model.n_bkwrd + model.n_both
        for i = 1:model.endo_nbr
            ws.a[i,model.i_bkwrd_b[j]] += ws.atmp[i,j]
        end
    end

    fill!(ws.b,0.0)
    b_cur = view(ws.b, :, model.i_fwrd_b)
    copy!(b_cur, ws.model_derivatives_1_fwrd)
    
    copy!(ws.c,view(first_order_ws.ghx,model.i_bkwrd_b,:))

    copy!(ws.ff, ws.model_derivatives_2_endo)
    A_mul_B!(ws.zz_fwrd, first_order_ws.ghx[model.i_fwrd_b,:], first_order_ws.ghx[model.i_bkwrd_b,:])
    copy!(ws.zz_current, first_order_ws.ghx)
    fill!(ws.zz_bkwrd, 0)
    println(size(ws.zz_bkwrd))
    for i=1:size(ws.zz, 2)
        ws.zz_bkwrd[i,i] = 1.0
    end
    
    order = 2
    display(ws.ff)
    display(ws.zz)
    a_mul_kron_b!(ws.d, ws.ff, ws.zz, order, ws.work)
    display(ws.d)
    eye_plus_at_kron_b_ws = EyePlusAtKronBWS(ws.a, ws.b, order, ws.c)

    vd = view(ws.d, :, 1:(model.n_bkwrd + model.n_both)^2)
    println("a")
    display(ws.a)
    println("b")
    display(ws.b)
    println("c")
    display(ws.c)
    println("d")
    display(ws.d)
    general_sylvester_solver!(ws.a, ws.b, ws.c, vd, order, eye_plus_at_kron_b_ws)
    println("d")
    display(ws.d)
end

end
