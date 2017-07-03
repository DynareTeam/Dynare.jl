module GeneralizedSylvester
###
# solving (I - s^T ⊗ s^T ⊗ ... \otimes s^T \otimes t)x = d
###

import ..LinAlg.QUT: QuasiUpperTriangular, A_mul_B!, I_plus_rA_ldiv_B!
import ..LinAlg.Kronecker: a_mul_b_kron_c!, kron_at_kron_b_mul_c!
import ..LinAlg.linsolve_algo: LinSolveWS, linsolve_core!
using  ..LinAlg.Schur: DgeesWS, dgees!

import Base.LinAlg.BLAS: gemm!
export EyePlusAtKronBWS, general_sylvester_solver!, real_eliminate!, solvi, transformation1

immutable EyePlusAtKronBWS
    n_s::Int64
    n_t::Int64
    vs_b::Matrix{Float64}
    vs_c::Matrix{Float64}
    s2::QuasiUpperTriangular{Float64,Matrix{Float64}}
    t2::QuasiUpperTriangular{Float64,Matrix{Float64}}
    work1::Vector{Float64}
    work2::Vector{Float64}
    work3::Vector{Float64}
    linsolve_ws::LinSolveWS
    dgees_ws_b::DgeesWS
    dgees_ws_c::DgeesWS
    function EyePlusAtKronBWS(a::AbstractArray,b::AbstractArray,order::Int64,c::AbstractArray)
        n_t, n_s = LinAlg.checksquare(a, b)
        mc, nc = size(c)
        if mc != n_t
            DimensionMismatch("a has dimensions ($n_t, $n_t) but c has dimensions ($mc, $nc)")
        elseif nc != nc^order
            DimensionMismatch("b has dimensions ($n_s, $n_s), order is $order, but c has dimensions  ($mc, $nc) and $n_s^$order = $(n_s^order)")
        end
        vs_b = Matrix{Float64}(n_t,n_t)
        vs_c = Matrix{Float64}(nc,nc)
        s2 = QuasiUpperTriangular(Matrix{Float64}(nc,nc))
        t2 = QuasiUpperTriangular(Matrix{Float64}(n_t,n_t))
        linsolve_ws = LinSolveWS(n_s)
        dgees_ws_b = DgeesWS(b)
        dgees_ws_c = DgeesWS(c)
        work1 = Vector{Float64}(n_t*nc^order)
        work2 = Vector{Float64}(n_t*nc^order)
        work3 = Vector{Float64}(n_t*nc^order)
        new(n_s,n_t,vs_b,vs_c,s2,t2,work1,work2,work3,linsolve_ws,dgees_ws_b,dgees_ws_c)
    end
end

"""
    function general_sylvester_solver!(a::AbstractMatrix,b::AbstractMatrix,c::AbstractMatrix,
                                       d::AbstractMatrix,order::Int64,ws::EyePlusAtKronBWS)
solves a*x + b*x*(c⊗c) = d
"""
function general_sylvester_solver!(a::AbstractMatrix,b::AbstractMatrix,c::AbstractMatrix,
                                   d::AbstractMatrix,order::Int64,ws::EyePlusAtKronBWS)
    linsolve_core!(ws.linsolve_ws,Ref{UInt8}('N'),a,b,d)
    dgees!(ws.dgees_ws_b,b)
    dgees!(ws.dgees_ws_c,c)
    t = QuasiUpperTriangular(b)
    A_mul_B!(ws.t2,t,t)
    s = QuasiUpperTriangular(c)
    A_mul_B!(ws.s2,s,s)
    result = reshape(ws.work1,size(t,2),size(s,1)^order)
    # do transpose in function
    a_mul_b_kron_c!(result,ws.dgees_ws_b.vs',d,ws.dgees_ws_c.vs,order)
    copy!(d,ws.work1)
    d_vec = vec(d)
    solve1!(1.0,order,t,ws.t2,s,ws.s2,d_vec,ws)
    # do transpose in function
    a_mul_b_kron_c!(result,ws.dgees_ws_b.vs,d,ws.dgees_ws_c.vs',order)
    copy!(d,result)
end


function solver!(t::QuasiUpperTriangular,s::QuasiUpperTriangular,d::AbstractVector,order::Int64,ws::EyePlusAtKronBWS)
    s2 = QuasiUpperTriangular(s*s)
    t2 = QuasiUpperTriangular(t*t)
    solve1!(1.0,order,t,t2,s,s2,d,ws)
    d
end

function solve1!(r,depth,t,t2,s,s2,d,ws)
    m = size(t,2)
    n = size(s,1)
    if depth == 0
        I_plus_rA_ldiv_B!(r,t,d)
    else
        nd = m*n^(depth-1)
        nd2 = 2*nd
        drange1 = 1:nd
        drange2 = 1:nd2
        i = 1
        while i <= n
            if i == n || s[i+1,i] == 0
                dv = view(d,drange1)
                solve1!(r*s[i,i],depth-1,t,t2,s,s2,dv,ws)
                if i < n
                    solvi_real_eliminate!(i,n,nd,drange1,depth-1,r,t,s,d,ws)
                end
                drange1 += nd
                drange2 += nd
                i += 1
            else
                dv = view(d,drange2)
                # s is transposed!
                solvii(r*s[i,i],r*s[i+1,i],r*s[i,i+1],depth-1,t,t2,s,s2,dv,ws)
                if i < n - 1
                    solvi_complex_eliminate!(i,n,nd,drange1,depth-1,r,t,s,d,ws)
                end
                drange1 += nd2
                drange2 += nd2
                i += 2
            end
        end
    end
end

function solvi_real_eliminate!(i::Int64,n::Int64,nd::Int64,drange::UnitRange{Int64},
                               depth::Int64,r::Float64,t::QuasiUpperTriangular,s::QuasiUpperTriangular,
                               d::AbstractVector,ws::EyePlusAtKronBWS)
    d_copy = view(ws.work1,1:length(d))
    copy!(d_copy,d)
    d1 = view(d_copy,drange)
    x = view(ws.work2,1:length(drange))
    kron_at_kron_b_mul_c!(x,s,depth,t,d1)
    @inbounds for j = i+1:n
        drange += nd
        dv = view(d,drange)
        m = r*s[i,j]
        @simd for k = 1:nd
            dv[k] -= m*x[k]
        end
    end 
end

function solvi_complex_eliminate!(i::Int64,n::Int64,nd::Int64,drange::UnitRange{Int64},
                                  depth::Int64,r::Float64,t::QuasiUpperTriangular,s::QuasiUpperTriangular,
                                  d::Vector{Float64},ws::EyePlusAtKronBWS)
    d_copy = view(ws.work1,1:length(d))
    copy!(d_copy,d)
    y1 =  view(ws.work2,1:length(drange))
    y2 =  view(ws.work3,1:length(drange))
    x1 = view(d_copy,drange)
    drange += nd
    x2 = view(d_copy,drange)
    kron_at_kron_b_mul_c!(y1,s,depth,t,x1)
    kron_at_kron_b_mul_c!(y2,s,depth,t,x2)
    @inbounds for j = i+2:n
        drange += nd
        dv = view(d,drange)
        # can probably be optimized
        m1 = r*s[i,j]
        m2 = r*s[i+1,j]
        @simd for k = 1:nd
            dv[k] -= m1*y1[k] + m2*y2[k]
        end
    end 
end

function solvii(alpha::Float64,beta1::Float64,beta2::Float64,depth::Int64,
                t::QuasiUpperTriangular,t2::QuasiUpperTriangular,s::QuasiUpperTriangular,
                s2::QuasiUpperTriangular,d::Vector{Float64},ws::EyePlusAtKronBWS)
    m = size(t,2)
    n = size(s,1)
    nd = m*n^depth
    transformation1(alpha,beta1,beta2,depth,t,s,d,ws)
    dv = view(d,1:nd)
    solviip(alpha,sqrt(-beta1*beta2),depth,t,t2,s,s2,dv,ws)
    dv = view(d,nd+1:2*nd)
    solviip(alpha,sqrt(-beta1*beta2),depth,t,t2,s,s2,dv,ws)
end

function transformation1(a::Float64,b1::Float64,b2::Float64,depth::Int64,
                         t::QuasiUpperTriangular,s::QuasiUpperTriangular,
                         d::Vector{Float64},ws::EyePlusAtKronBWS)
    m = size(t,2)
    n = size(s,1)
    nd = m*n^depth
    d_orig = view(ws.work1,1:length(d))
    copy!(d_orig,d)
    drange = 1:nd
    d1 = view(d_orig,drange)
    d2 = view(d_orig,drange+nd)
    work = view(ws.work2,drange)
    kron_at_kron_b_mul_c!(s,depth,t,d1,work)
    kron_at_kron_b_mul_c!(s,depth,t,d2,work)
    @inbounds @simd for i = drange
        d[i] += a*d1[i] - b1*d2[i]
        d[i+nd] += -b2*d1[i] + a*d2[i] 
    end
end

diag_zero_sq = 1e-30

function solviip(alpha::Float64,beta::Float64,depth::Int64,t::QuasiUpperTriangular,t2::QuasiUpperTriangular,
                 s::QuasiUpperTriangular,s2::QuasiUpperTriangular,d::Vector{Float64},ws::EyePlusAtKronBWS)
    m = size(t,2)
    n = size(s,1)
    if beta*beta < diag_zero_sq
        solve1!(alpha,depth,t,t2,s,s2,d,ws)
        solve1!(alpha,depth,t,t2,s,s2,d,ws)
        return
    end

    if depth == 0
        I_plus_rA_plus_sB_ldiv_C!(2*alpha,alpha*alpha+beta*beta,t,t2,d)
    else
        nd = m*n^(depth-1)
        nd2 = 2*nd
        drange1 = 1:nd
        drange2 = 1:nd2
        i = 1
        while i <= n
            # s is transposed
            if i == n || s[i+1,i] == 0
                dv = view(d,drange1)
                if s[i,i]*s[i,i]*(alpha*alpha+beta*beta) > diag_zero_sq
                    solviip(s[i,i]*alpha,s[i,i]*beta,depth-1,t,t2,s,s2,dv,ws)
                end
                if i < n
                    solviip_real_eliminate!(i,n,nd,drange1,depth-1,alpha,beta,t,t2,s,s2,d,ws)
                end
                drange1 += nd
                drange2 += nd
                i += 1
            else
                dv = view(d,drange2)
                # s transposed !
                solviip2(alpha,beta,s[i,i],s[i+1,i],s[i,i+1],depth,t,t2,s,s2,dv,ws)
                if i < n - 1
                    solviip_complex_eliminate!(i,n,nd,drange1,depth-1,alpha,beta,t,t2,s,s2,d,ws)
                end
                drange1 += nd2
                drange2 += nd2
                i += 2
            end
        end
    end
end

function solviip_real_eliminate!(i::Int64,n::Int64,nd::Int64,drange::UnitRange{Int64},
                                 depth::Int64,alpha::Float64,beta::Float64,t::QuasiUpperTriangular,
                                 t2::QuasiUpperTriangular,s::QuasiUpperTriangular,
                                 s2::QuasiUpperTriangular,d::Vector{Float64},ws::EyePlusAtKronBWS)
    y1 = view(ws.work1,drange)
    y2 = view(ws.work1,drange+nd)
    copy!(y1,d[drange])
    copy!(y2,d[drange])
    work = view(ws.work2,1:length(drange))
    kron_at_kron_b_mul_c!(s,depth,t,y1,work)
    kron_at_kron_b_mul_c!(s2,depth,t2,y2,work)
    @inbounds for j = i+1:n
        drange += nd
        dv = view(d,drange)
        m1 = 2*alpha*s[i,j]
        m2 = (alpha*alpha+beta*beta)*s2[i,j]
        @simd for k = 1:nd
            dv[k] -= m1*y1[k] + m2*y2[k]
        end
    end 
end

function solviip2(alpha::Float64,beta::Float64,gamma::Float64,delta1::Float64,delta2::Float64,
                  depth::Int64,t::QuasiUpperTriangular,t2::QuasiUpperTriangular,
                  s::QuasiUpperTriangular,s2::QuasiUpperTriangular,d::Vector{Float64},ws::EyePlusAtKronBWS)
    m = size(t,2)
    n = size(s,1)
    aspds = alpha*alpha + beta*beta
    gspds = gamma*gamma - delta1*delta2
    nd = m*n^(depth-1)
    dv1 = view(d,1:nd)
    dv2 = view(d,nd+1:2*nd)
    if aspds*gspds > diag_zero_sq
        transform2(alpha, beta, gamma, -delta1, -delta2, nd, depth, t, t2, s, s2, d,ws)
        delta = sqrt(-delta1*delta2)
	a1 = alpha*gamma - beta*delta
	b1 = alpha*delta + gamma*beta
	a2 = alpha*gamma + beta*delta
	b2 = alpha*delta - gamma*beta
	solviip(a2, b2, depth-1, t, t2, s, s2, dv1,ws);
	solviip(a1, b1, depth-1, t, t2, s, s2, dv1,ws);
	solviip(a2, b2, depth-1, t, t2, s, s2, dv2,ws);
	solviip(a1, b1, depth-1, t, t2, s, s2, dv2,ws);
    end
end

function transform2(alpha::Float64, beta::Float64, gamma::Float64, delta1::Float64, delta2::Float64,
                    nd::Float64, depth::Int64, t::QuasiUpperTriangular, t2::QuasiUpperTriangular,
                    s::QuasiUpperTriangular, s2::QuasiUpperTriangular, d::Vector{Float64}, ws::EyePlusAtKronBWS)
    drange = 1:nd
    d1 = view(d,drange)
    d2 = view(d,drange+nd)
    d1tmp = view(ws.work1,drange)
    d2tmp = view(ws.work1,drange+nd)
    copy!(d1tmp,d1)
    copy!(d2tmp,d2)
    x1 = view(ws.work2,drange)
    x2 = view(ws.work2,drange+nd)
    copy!(x1,d1);
    copy!(x2,d2);
    work = view(ws.work3,1:length(drange))
    kron_at_kron_b_mul_c!(s,depth-1,t,d1tmp,work)
    kron_at_kron_b_mul_c!(s,depth-1,t,d2tmp,work)

    m1 = 2*alpha*gamma
    m2 = 2*alpha*delta1
    @inbounds @simd for i = 1:nd
        d1[i] += m1*d1tmp[i] + m2*d2tmp[i]
        d2[i] += 2*alpha*(delta2*d1tmp[i] + gamma*d2tmp[i])
    end
    
    copy!(d1tmp,x1); # restore to d1
    copy!(d2tmp,x2); # restore to d2
    kron_at_kron_b_mul_c!(s2,depth-1,t2,d1tmp,work)
    kron_at_kron_b_mul_c!(s2,depth-1,t2,d2tmp,work)

    aspds = alpha*alpha + beta*beta;
    gspds = gamma*gamma + delta1*delta2;
    m1 = aspds*gspds
    m2 = 2*aspds*gamma*delta1
    m3 = 2*aspds*gamma*delta2
    @inbounds @simd for i = 1:nd
        d1[i] += m1*d1tmp[i] + m2*d2tmp[i]
        d2[i] += m3*d1tmp[i] + m1*d2tmp[i]
    end
end

"""
    solviip_complex_eliminate!(i,n,nd,drange,depth,alpha,beta,t,t2,s,s2,d)

perfoms elimination after solving for a complex diagonal block of size 2*n^depth

d n^(depth+2) x 1

The solution is stored in d[drange; drange + nd]

The function updates d[i*nd+1:n*nd]
"""
function solviip_complex_eliminate!(i::Int64,n::Int64,nd::Int64,drange::UnitRange{Int64},depth::Int64,
                                    alpha::Float64,beta::Float64,t::QuasiUpperTriangular,t2::QuasiUpperTriangular,
                                    s::QuasiUpperTriangular,s2::QuasiUpperTriangular,d::Vector{Float64},ws::EyePlusAtKronBWS)
    y11 = view(ws.work1,drange)
    y12 = view(ws.work1,drange+nd)
    copy!(y11,d[drange])
    copy!(y12,d[drange])
    drange += nd
    y21 = view(ws.work2,drange)
    y22 = view(ws.work2,drange+nd)
    copy!(y21,d[drange])
    copy!(y22,d[drange])
    work = view(ws.work3,drange)

    kron_at_kron_b_mul_c!(s,depth,t,y11,work)
    kron_at_kron_b_mul_c!(s2,depth,t2,y12,work)
    kron_at_kron_b_mul_c!(s,depth,t,y21,work)
    kron_at_kron_b_mul_c!(s2,depth,t2,y22,work)

    alpha2beta2 = alpha*alpha + beta*beta
    @inbounds for j = i+2:n
        drange += nd
        dv = view(d,drange)
        m1 = 2*alpha*s[i,j]
        m2 = alpha2beta2*s2[i,j]
        m3 = 2*alpha*s[i+1,j]
        m4 = alpha2beta2*s2[i+1,j]
        @simd for k = 1:nd
            dv[k] -= m1*y11[k] + m2*y12[k] + m3*y21[k] + m4*y22[k]
        end
    end 
end

end
