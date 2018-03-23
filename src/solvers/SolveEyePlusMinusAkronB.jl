module SolveEyePlusMinusAkronB
###
# solving a x + b x (c ⊗ c ⊗ ... ⊗ c) = d
# using (I - s^T ⊗ s^T ⊗ ... \otimes s^T \otimes t)x = d'
###

using ...DynLinAlg.QUT
using ...DynLinAlg.LinSolveAlgo
using ...DynLinAlg.SchurAlgo
using ...DynLinAlg.KroneckerUtils
using Base.Test
import Base.LinAlg.BLAS: gemm!
export EyePlusAtKronBWS, generalized_sylvester_solver!, real_eliminate!, solvi, transformation1

immutable EyePlusAtKronBWS
    ma::Int64
    mb::Int64
    b1::Matrix{Float64}
    c1::Matrix{Float64}
    vs_b::Matrix{Float64}
    vs_c::Matrix{Float64}
    s2::QuasiUpperTriangular{Float64,Matrix{Float64}}
    t2::QuasiUpperTriangular{Float64,Matrix{Float64}}
    work1::Vector{Float64}
    work2::Vector{Float64}
    work3::Vector{Float64}
    result::Matrix{Float64}
    linsolve_ws::LinSolveWS
    dgees_ws_b::DgeesWS
    dgees_ws_c::DgeesWS
    function EyePlusAtKronBWS(ma::Int64, mb::Int64, mc::Int64, order::Int64)
        if mb != ma
            DimensionMismatch("a has $ma rows but b has $mb rows")
        end
        b1 = Matrix{Float64}(mb,mb)
        c1 = Matrix{Float64}(mc,mc)
        vs_b = Matrix{Float64}(mb,mb)
        vs_c = Matrix{Float64}(mc,mc)
        s2 = QuasiUpperTriangular(Matrix{Float64}(mc,mc))
        t2 = QuasiUpperTriangular(Matrix{Float64}(mb,mb))
        linsolve_ws = LinSolveWS(ma)
        dgees_ws_b = DgeesWS(mb)
        dgees_ws_c = DgeesWS(mc)
        work1 = Vector{Float64}(ma*mc^order)
        work2 = Vector{Float64}(ma*mc^order)
        work3 = Vector{Float64}(ma*mc^order)
        result = Matrix{Float64}(ma,mc^order)
        new(ma, mb, b1, c1, vs_b, vs_c, s2, t2, work1, work2, work3, result, linsolve_ws, dgees_ws_b, dgees_ws_c)
    end
end

function generalized_sylvester_solver!(a::AbstractMatrix,b::AbstractMatrix,c::AbstractMatrix,
                                   d::AbstractVector,order::Int64,ws::EyePlusAtKronBWS)
    copy!(ws.b1,b)
    copy!(ws.c1,c)
    d1 = reshape(d,size(b,2),size(c,1)^order)
    linsolve_core!(ws.linsolve_ws,Ref{UInt8}('N'),a,ws.b1,d1)
    dgees!(ws.dgees_ws_b,ws.b1)
    dgees!(ws.dgees_ws_c,ws.c1)
    t = QuasiUpperTriangular(ws.b1)
    A_mul_B!(ws.t2,t,t)
    s = QuasiUpperTriangular(ws.c1)
    A_mul_B!(ws.s2,s,s)
    at_mul_b_kron_c!(ws.result, ws.dgees_ws_b.vs, d1, ws.dgees_ws_c.vs, order, ws.work2, ws.work3)
    copy!(d,ws.result)
    solve1!(1.0, order, t, ws.t2, s, ws.s2, d, ws)
    a_mul_b_kron_ct!(ws.result, ws.dgees_ws_b.vs, d1, ws.dgees_ws_c.vs, order, ws.work2, ws.work3)
    copy!(d,ws.result)
end


function solver!(t::QuasiUpperTriangular,s::QuasiUpperTriangular,d::AbstractVector,order::Int64,ws::EyePlusAtKronBWS)
    s2 = QuasiUpperTriangular(s*s)
    t2 = QuasiUpperTriangular(t*t)
    solve1!(1.0,order,t,t2,s,s2,d,ws)
    d
end

function solve1!(r::Float64, depth::Int64, t::AbstractArray{Float64,2}, t2::AbstractArray{Float64,2}, s::AbstractArray{Float64,2}, s2::AbstractArray{Float64,2}, d::AbstractVector{Float64}, ws::EyePlusAtKronBWS)
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
    work1 = ws.work1
    work2 = ws.work2
    work3 = ws.work3
    kron_at_kron_b_mul_c!(work1,1,s,depth,t,d,drange[1],work2,work3,1)
    k1 = drange[1] + nd
    @inbounds for j = i+1:n
        m = r*s[i,j]
        @simd for k2 = 1:nd
            d[k1] -= m*work1[k2]
            k1 += 1
        end
    end
end

function solvi_complex_eliminate!(i::Int64,n::Int64,nd::Int64,drange::UnitRange{Int64},
                                  depth::Int64,r::Float64,t::QuasiUpperTriangular,s::QuasiUpperTriangular,
                                  d::AbstractVector,ws::EyePlusAtKronBWS)
    work1 = ws.work1
    work2 = ws.work2
    work3 = ws.work3
    kron_at_kron_b_mul_c!(work1, 1, s, depth, t, d, drange[1], work2, work3, 1)
    drange += nd
    kron_at_kron_b_mul_c!(work1, nd+1, s, depth, t, d, drange[1], work2, work3, 1)
    k1 = drange[1] + nd
    @inbounds for j = i + 2 : n
        m1 = r*s[i,j]
        m2 = r*s[i+1,j]
        @simd for k2 = 1:nd
            d[k1] -= m1*work1[k2] + m2*work1[k2 + nd]
            k1 += 1
        end
    end 
end

function solvii(alpha::Float64,beta1::Float64,beta2::Float64,depth::Int64,
                t::QuasiUpperTriangular,t2::QuasiUpperTriangular,s::QuasiUpperTriangular,
                s2::QuasiUpperTriangular,d::AbstractVector,ws::EyePlusAtKronBWS)
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
                         d::AbstractVector,ws::EyePlusAtKronBWS)
    m = size(t,2)
    n = size(s,1)
    nd = m*n^depth
    copy!(ws.work3,d)
    drange = 1:nd
    d1 = view(ws.work3,drange)
    d2 = view(ws.work3,drange+nd)
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
                 s::QuasiUpperTriangular,s2::QuasiUpperTriangular,d::AbstractVector,ws::EyePlusAtKronBWS)
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
                                 s2::QuasiUpperTriangular,d::AbstractVector,ws::EyePlusAtKronBWS)
    y1 = view(ws.work1,drange)
    y2 = view(ws.work1,drange+nd)
    copy!(y1,1,d,drange[1],nd)
    copy!(y2,1,d,drange[1],nd)
    work = view(ws.work2,1:length(drange))
    kron_at_kron_b_mul_c!(s,depth,t,y1,work)
    kron_at_kron_b_mul_c!(s2,depth,t2,y2,work)
    k1 = drange[1] + nd
    @inbounds for j = i+1:n
        m1 = 2*alpha*s[i,j]
        m2 = (alpha*alpha+beta*beta)*s2[i,j]
        @simd for k2 = 1:nd
            d[k1] -= m1*y1[k2] + m2*y2[k2]
            k1 += 1
        end
    end 
end

function solviip2(alpha::Float64,beta::Float64,gamma::Float64,delta1::Float64,delta2::Float64,
                  depth::Int64,t::QuasiUpperTriangular,t2::QuasiUpperTriangular,
                  s::QuasiUpperTriangular,s2::QuasiUpperTriangular,d::AbstractVector,ws::EyePlusAtKronBWS)
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
	solviip(a2, b2, depth-1, t, t2, s, s2, dv1, ws);
	solviip(a1, b1, depth-1, t, t2, s, s2, dv1, ws);
        solviip(a2, b2, depth-1, t, t2, s, s2, dv2, ws);
        solviip(a1, b1, depth-1, t, t2, s, s2, dv2, ws);
    end
end

function transform2(alpha::Float64, beta::Float64, gamma::Float64, delta1::Float64, delta2::Float64,
                    nd::Int64, depth::Int64, t::QuasiUpperTriangular, t2::QuasiUpperTriangular,
                    s::QuasiUpperTriangular, s2::QuasiUpperTriangular, d::AbstractVector, ws::EyePlusAtKronBWS)
    d1 = ws.work1
    kron_at_kron_b_mul_c!(d1,1,s,depth-1,t,d,1,ws.work2,ws.work3,1)
    kron_at_kron_b_mul_c!(d1,nd+1,s,depth-1,t,d,nd+1,ws.work2,ws.work3,1)

    m1 = 2*alpha*gamma
    m2 = 2*alpha*delta1
    @inbounds @simd for i = 1:nd
        dtmp = d1[i]
        d1[i] = d[i] + m1*dtmp + m2*d1[i+nd]
        d1[i+nd] = d[i+nd] + 2*alpha*(delta2*dtmp + gamma*d1[i+nd])
    end

    kron_at_kron_b_mul_c!(d,1,s2,depth-1,t2,d,1,ws.work2,ws.work3,1)
    kron_at_kron_b_mul_c!(d,nd+1,s2,depth-1,t2,d,nd+1,ws.work2,ws.work3,1)
    
    aspds = alpha*alpha + beta*beta;
    gspds = gamma*gamma + delta1*delta2;
    m1 = aspds*gspds
    m2 = 2*aspds*gamma*delta1
    m3 = 2*aspds*gamma*delta2
    @inbounds @simd for i = 1:nd
        dtmp = d[i]
        d[i] = d1[i] + m1*dtmp + m2*d[i+nd]
        d[i+nd] = d1[i+nd] + m3*dtmp + m1*d[i+nd]
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
                                    s::QuasiUpperTriangular,s2::QuasiUpperTriangular,d::AbstractVector,ws::EyePlusAtKronBWS)
    y11 = view(ws.work1,drange)
    y12 = view(ws.work1,drange+nd)
    copy!(y11,1,d,drange[1],nd)
    copy!(y12,1,d,drange[1],nd)
    drange += nd
    y21 = view(ws.work2,drange)
    y22 = view(ws.work2,drange+nd)
    copy!(y21,1,d,drange[1],nd)
    copy!(y22,1,d,drange[1],nd)
    work = view(ws.work3,drange)

    kron_at_kron_b_mul_c!(s,depth,t,y11,work)
    kron_at_kron_b_mul_c!(s2,depth,t2,y12,work)
    kron_at_kron_b_mul_c!(s,depth,t,y21,work)
    kron_at_kron_b_mul_c!(s2,depth,t2,y22,work)

    alpha2beta2 = alpha*alpha + beta*beta
    k1 = drange[1] + nd
    @inbounds for j = i+2:n
        m1 = 2*alpha*s[i,j]
        m2 = alpha2beta2*s2[i,j]
        m3 = 2*alpha*s[i+1,j]
        m4 = alpha2beta2*s2[i+1,j]
        @simd for k2 = 1:nd
            d[k1] -= m1*y11[k2] + m2*y12[k2] + m3*y21[k2] + m4*y22[k2]
            k1 += 1
        end
    end 
end

end
