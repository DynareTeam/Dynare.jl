module CyclicReduction

import ...DynLinAlg.LinSolveAlgo: LinSolveWS, linsolve_core!

import Base.LinAlg.BLAS: scal!, gemm!

export CyclicReductionWS, cyclic_reduction!, cyclic_reduction_check

type CyclicReductionWS
    linsolve_ws::LinSolveWS
    ahat1::Matrix{Float64}
    a1copy::Matrix{Float64}
    m::Matrix{Float64,}
    m00::SubArray{Float64}
    m02::SubArray{Float64}
    m20::SubArray{Float64}
    m22::SubArray{Float64}
    m1::Matrix{Float64}
    m1_a0::SubArray{Float64}
    m1_a2::SubArray{Float64}
    m2::Matrix{Float64}
    m2_a0::SubArray{Float64}
    m2_a2::SubArray{Float64}
    info::Int64

    function CyclicReductionWS(n)
        linsolve_ws = LinSolveWS(n)
        ahat1 = Matrix{Float64}(n,n)
	a1copy = Matrix{Float64}(n,n)
        m = Matrix{Float64}(2*n,2*n)
        m00 = view(m,1:n,1:n)
        m02 = view(m,1:n,n+(1:n))
        m20 = view(m,n+(1:n),1:n)
        m22 = view(m,n+(1:n),n+(1:n))
        m1 = Matrix{Float64}(n,2*n)
        m1_a0 = view(m1,1:n,1:n)
        m1_a2 = view(m1,1:n,n+(1:n))
        m2 = Matrix{Float64}(2*n,n)
        m2_a0 = view(m2,1:n,1:n)
        m2_a2 = view(m2,n+(1:n),1:n)
        info = 0
        new(linsolve_ws,ahat1,a1copy,m, m00, m02, m20, m22, m1, m1_a0, m1_a2, m2, m2_a0, m2_a2, info) 
    end
end

"""
    cyclic_reduction!(x::Array{Float64},a0::Array{Float64},a1::Array{Float64},a2::Array{Float64},ws::CyclicReductionWS, cvg_tol::Float64, max_it::Int64)

Solve the quadratic matrix equation a0 + a1*x + a2*x*x = 0, using the cyclic reduction method from Bini et al. (???).

The solution is returned in matrix x. In case of nonconvergency, x is set to NaN and an error code is returned in ws.info
* info = 0: return OK
* info = 1: no stable solution (????)
* info = 2: multiple stable solutions (????)

# Example
```meta
DocTestSetup = quote
     using CyclicReduction
     n = 3
     ws = CyclicReductionWS(n)
     a0 = [0.5 0 0; 0 0.5 0; 0 0 0];
     a1 = eye(n)
     a2 = [0 0 0; 0 0 0; 0 0 0.8]
     x = zeros(n,n)
end
```

```jldoctest
julia> display(names(CyclicReduction))
```

```jldoctest
julia> cyclic_reduction!(x,a0,a1,a2,ws,1e-8,50)
```
"""
function cyclic_reduction!(x::Array{Float64},a0::Array{Float64},a1::Array{Float64},a2::Array{Float64},ws::CyclicReductionWS, cvg_tol::Float64, max_it::Int64)
    copy!(x,a0)
    copy!(ws.ahat1,1,a1,1,length(a1))
    @inbounds copy!(ws.m1_a0, a0)
    @inbounds copy!(ws.m1_a2, a2)
    @inbounds copy!(ws.m2_a0, a0)
    @inbounds copy!(ws.m2_a2, a2)
    it = 0
    @inbounds while it < max_it
        #        ws.m = [a0; a2]*(a1\[a0 a2])
	copy!(ws.a1copy,a1)
        linsolve_core!(ws.linsolve_ws,Ref{UInt8}('N'),ws.a1copy,ws.m1)
        gemm!('N','N',-1.0,ws.m2,ws.m1,0.0,ws.m)
        @simd for i in eachindex(a1)
            a1[i] += ws.m02[i] + ws.m20[i]
        end
        copy!(ws.m1_a0, ws.m00)
        copy!(ws.m1_a2, ws.m22)
        copy!(ws.m2_a0, ws.m00)
        copy!(ws.m2_a2, ws.m22)
        if any(isinf.(ws.m))
            if norm(ws.m1_a0) < cvg_tol
                ws.info = 2
            else
                ws.info = 1
            end
            fill!(x,NaN)
            return
        end
        ws.ahat1 += ws.m20
        crit = norm(ws.m1_a0,1)
        if crit < cvg_tol
	    # keep iterating until condition on a2 is met
            if norm(ws.m1_a2,1) < cvg_tol
                break
            end
        end
        it += 1
    end
    if it == max_it
        if norm(ws.m1_a0) < cvg_tol
            ws.info = 2
        else
            ws.info = 1
        end
        fill!(x,NaN)
        return
    else
        linsolve_core!(ws.linsolve_ws,Ref{UInt8}('N'),ws.ahat1,x)
        @inbounds scal!(length(x),-1.0,x,1)
        ws.info = 0
    end        
end

function cyclic_reduction_check(x::Array{Float64,2},a0::Array{Float64,2}, a1::Array{Float64,2}, a2::Array{Float64,2},cvg_tol::Float64)
    res = a0 + a1*x + a2*x*x
    if (sum(sum(abs.(res))) > cvg_tol)
        print("the norm of the residuals, ", res, ", compared to the tolerance criterion ",cvg_tol)
    end
    nothing
end
     
end
