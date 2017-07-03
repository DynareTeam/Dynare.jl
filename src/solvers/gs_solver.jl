include("exceptions.jl")

module gs_solver

import ..LinAlg.Schur: DggesWS, dgges!
import ..LinAlg.linsolve_algo: LinSolveWS, linsolve_core!

import Base.LinAlg.BLAS: scal!, gemm!

export GsSolverWS, gs_solver!

type GsSolverWS
    dgges_ws::DggesWS
    linsolve_ws::LinSolveWS
    D11::SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}
    E11::SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}
    Z11::SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}
    Z12::SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}
    Z21::SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}
    Z22::SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}
    tmp1::Array{Float64,2}
    tmp2::Array{Float64,2}
    tmp3::Array{Float64,2}
    g1::Array{Float64,2}
    g2::Array{Float64,2}
#    gx_f::SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},Colon},false}
#    hx::Array{Float64,2}
#    ghx::Array{Float64,2}
#    ghx_bb::SubArray{Float64,2,Array{Float64,2},Tuple{Array{Int64,1},Colon},false}
#    ghx_f::SubArray{Float64,2,Array{Float64,2},Tuple{Array{Int64,1},Colon},false}

    function GsSolverWS(d,e,n1)
        dgges_ws = DggesWS(Ref{UInt8}('N'),Ref{UInt8}('V'),e,d)
        n = size(d,1)
        n2 = n - n1
        linsolve_ws = LinSolveWS(n2)
        D11 = view(d,1:n1,1:n1)
        E11 = view(e,1:n1,1:n1)
        Z11 = view(dgges_ws.vsr,1:n1,1:n1)
        Z12 = view(dgges_ws.vsr,1:n1,n1+1:n)
        Z21 = view(dgges_ws.vsr,n1+1:n,1:n1)
        Z22 = view(dgges_ws.vsr,n1+1:n,n1+1:n)
        tmp1 = Array(Float64,n2,n1)
        tmp2 = Array(Float64,n1,n1)
        tmp3 = Array(Float64,n1,n1)
        g1 = Array(Float64,n1,n1)
        g2 = Array(Float64,n2,n1)
#        gx_f = view(gx,1:n2,:)
#        ghx = Array(Float64,model.endo_nbr,n1)
#        ghx_bb = view(ghx,i1,:)
#        ghx_f = view(ghx,i2,:)
        new(dgges_ws,linsolve_ws,D11,E11,Z11,Z12,Z21,Z22,tmp1,tmp2,tmp3,g1,g2)
    end
end

"""
    gs_solver!(ws::GsSolverWS,d::Array{Float64,2},e::Array{Float64,2},n1::Int64,qz_criterium)

finds the unique stable solution for the following system:

```
d \left[\begin{array}{c}I\\g_2\end{array}\right]g_1 = e \left[\begin{array}{c}I\\g_2\end{array}\right]
```
"""
function gs_solver!(ws::GsSolverWS,d::Array{Float64,2},e::Array{Float64,2},n1::Int64,qz_criterium::Float64)

    dgges!(ws.dgges_ws,e,d)
    nstable = ws.dgges_ws.sdim[]
    
    if nstable < n1
        throw(UnstableSystemException)
    elseif nstable > n1
        throw(UndeterminateSystemExcpetion)
    end
    ws.g2 = ws.Z12'
#    ws.tmp2 = ws.Z22'
#    gx = -(ws.Z12/ws.Z22)'
    linsolve_core!(ws.linsolve_ws,Ref{UInt8}('T'),ws.Z22,ws.g2)
    scal!(length(ws.g2),-1.0,ws.g2,1)
    ws.tmp2 = ws.Z11'
#    hx1 = ws.dgges_ws.vsr[1:nstable,1:nstable]/D[1:nstable, 1:nstable]
    ws.D11 = view(d,1:nstable,1:nstable)
    linsolve_core!(ws.linsolve_ws,Ref{UInt8}('T'),ws.D11,ws.tmp2)
#    hx2 = E[1:nstable,1:nstable]/ws.dgges_ws.vsr[1:nstable,1:nstable]
    ws.E11 = view(e,1:nstable,1:nstable)
    ws.tmp3 = ws.E11'
    linsolve_core!(ws.linsolve_ws,Ref{UInt8}('T'),ws.Z11,ws.tmp3)
    #hx = ws.tmp2*ws.tmp3
    gemm!('T','T',1.0,ws.tmp2,ws.tmp3,0.0,ws.g1)
#    if false
#        for i = 1:ns
#            for j = 1:model.n_bkwrd
#                ws.ghx[model.i_bkwrd[j],i] = ws.hx[j,i]
#            end
#            for j = 1:model.n_fwrd
#                ws.ghx[model.i_fwrd[j],i] = ws.g2[j,i]
#            end
#            for j = 1:model.n_both
#                ws.ghx[model.i_both[j],i] = ws.hx[model.n_bkwrd+j,i]
#            end
#        end
#    else
#        ws.ghx[model.i_bkwrd,:] = ws.hx[1:model.n_bkwrd,:]
#        ws.ghx[model.i_fwrd,:] = ws.g2[1:model.n_fwrd,:]
#        ws.ghx[model.i_both,:] = ws.g2[model.n_fwrd+(1:model.n_both),:]
#    end
end

end
