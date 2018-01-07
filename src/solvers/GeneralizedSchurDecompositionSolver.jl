include("exceptions.jl")

module GeneralizedSchurDecompositionSolver

import ...DynLinAlg.SchurAlgo: DggesWS, dgges!
import ...DynLinAlg.LinSolveAlgo: LinSolveWS, linsolve_core!

import Base.LinAlg.BLAS: scal!, gemm!

export GsSolverWS, gs_solver!

type GsSolverWS
    dgges_ws::DggesWS
    linsolve_ws::LinSolveWS
    D11::SubArray{Float64}
    E11::SubArray{Float64}
    Z11::SubArray{Float64}
    Z12::SubArray{Float64}
    Z21::SubArray{Float64}
    Z22::SubArray{Float64}
    tmp1::Matrix{Float64}
    tmp2::Matrix{Float64}
    tmp3::Matrix{Float64}
    g1::Matrix{Float64}
    g2::Matrix{Float64}
    eigval::Vector{Complex64}
    
    function GsSolverWS(d,e,n1)
        dgges_ws = DggesWS(Ref{UInt8}('N'), Ref{UInt8}('N'), Ref{UInt8}('N'), e, d)
        n = size(d,1)
        n2 = n - n1
        linsolve_ws = LinSolveWS(n2)
        D11 = view(d,1:n1,1:n1)
        E11 = view(e,1:n1,1:n1)
        vsr = Matrix{Float64}(n,n)
        Z11 = view(vsr,1:n1,1:n1)
        Z12 = view(vsr,1:n1,n1+1:n)
        Z21 = view(vsr,n1+1:n,1:n1)
        Z22 = view(vsr,n1+1:n,n1+1:n)
        tmp1 = Matrix{Float64}(n2,n1)
        tmp2 = Matrix{Float64}(n1,n1)
        tmp3 = Matrix{Float64}(n1,n1)
        g1 = Matrix{Float64}(n1,n1)
        g2 = Matrix{Float64}(n2,n1)
        eigval = Vector{Complex64}(n)
        new(dgges_ws,linsolve_ws,D11,E11,Z11,Z12,Z21,Z22,tmp1,tmp2,tmp3,g1,g2, eigval)
    end
end

"""
    gs_solver!(ws::GsSolverWS,d::Matrix{Float64},e::Matrix{Float64},n1::Int64,qz_criterium)

finds the unique stable solution for the following system:

```
d \left[\begin{array}{c}I\\g_2\end{array}\right]g_1 = e \left[\begin{array}{c}I\\g_2\end{array}\right]
```
"""
function gs_solver!(ws::GsSolverWS,d::Matrix{Float64},e::Matrix{Float64},n1::Int64,qz_criterium::Float64)

    dgges!('N', 'V', e, d, [0.0], ws.vsr, ws.eigval, ws.dgges_ws)
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
