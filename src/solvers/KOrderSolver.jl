module KOrderSolver

using Base.Test
using ...DynLinAlg.KroneckerUtils
using ...DynLinAlg.LinSolveAlgo
import Base.LinAlg.BLAS: gemm!
import ...Solvers.SolveEyePlusMinusAkronB: EyePlusAtKronBWS, generalized_sylvester_solver!
import ...FaaDiBruno: partial_faa_di_bruno!, FaaDiBrunoWs
export make_gg!, make_hh!, k_order_solution!, KOrderWs

struct KOrderWs
    nvar::Integer
    nfwrd::Integer
    nstate::Integer
    ncur::Integer
    nshock::Integer
    fwrd_index::Array{Int64}
    state_index::Array{Int64}
    cur_index::Array{Int64}
    state_range::Range
    gfwrd::Vector{Matrix{Float64}}
    gg::Vector{Matrix{Float64}}
    hh::Vector{Matrix{Float64}}
    rhs::Matrix{Float64}
    rhs1::Matrix{Float64}
    gykf::Matrix{Float64}
    gs_su::Matrix{Float64}
    a::Matrix{Float64}
    b::Matrix{Float64}
    work1::Vector{Float64}
    work2::Vector{Float64}
    faa_di_bruno_ws_1::FaaDiBrunoWs
    faa_di_bruno_ws_2::FaaDiBrunoWs
    linsolve_ws_1::LinSolveWS
    gs_ws::EyePlusAtKronBWS
    function KOrderWs(nvar,nfwrd,nstate,ncur,nshock,fwrd_index,state_index,cur_index,state_range,order)
        gfwrd = [zeros(nfwrd,(nstate+nshock+1)^i) for i = 1:order]
        gg = [zeros(nstate+nshock+1,(nstate+2*nshock+1)^i) for i = 1:order]
        hh = [zeros(nfwrd+nvar+nstate+nshock,(nstate+2*nshock+1)^i) for i = 1:order]
        faa_di_bruno_ws_1 = FaaDiBrunoWs(nfwrd, nstate + 2*nshock + 1, order)
        faa_di_bruno_ws_2 = FaaDiBrunoWs(nvar, nfwrd + ncur + nstate + nshock, order)
        linsolve_ws_1 = LinSolveWS(nvar)
        rhs = zeros(nvar,(nstate+2*nshock+1)^order)
        rhs1 = zeros(nvar, max(nvar^order,nshock*(nstate+nshock)^(order-1)))
        gykf = zeros(nfwrd,nstate^order)
        gs_su = Array{Float64}(nstate,nstate+nshock)
        a = zeros(nvar,nvar)
        b = zeros(nvar,nvar)
        work1 = zeros(nvar*(nstate + nshock + 1)^order)
        work2 = similar(work1)
        gs_ws = EyePlusAtKronBWS(nvar,nvar,nstate,order)
        new(nvar,nfwrd,nstate,ncur,nshock,fwrd_index,state_index,cur_index,state_range,gfwrd,gg,hh,
            rhs,rhs1,gykf,gs_su,a,b,work1,work2,faa_di_bruno_ws_1,faa_di_bruno_ws_2,linsolve_ws_1, gs_ws)
    end
end

"""
    function make_gg!(gg,g,order,ws)

assembles the derivatives of function
gg(y,u,σ,ϵ) = [g_state(y,u,σ); ϵ; σ] at order 'order' 
with respect to [y, u, σ, ϵ]
"""  
function make_gg!(gg,g,order,ws)
    ngg1 = ws.nstate + 2*ws.nshock + 1
    mgg1 = ws.nstate + ws.nshock + 1
    @assert size(gg[order]) == (mgg1, ngg1^order)
    @assert size(g[order],2) == (ws.nstate + ws.nshock + 1)^order
    @assert ws.state_range.stop <= size(g[order],1)
    if order == 1
        v2 = view(g[1],ws.state_index,:)
        copy!(gg[1],v2)
        for i = 1:ws.nshock
            gg[1][ws.nstate + i, ws.nstate + ws.nshock + 1 + i] = 1.0
        end
        gg[1][end, ws.nstate + ws.nshock + 1] = 1 
    else
        n = ws.nstate + ws.nshock + 1
        i1 = CartesianIndex(1,(repmat([1], order - 1))...)
        i2 = CartesianIndex(1,(repmat([n], order - 1))...)
        i3 = ((repmat([ngg1], order))...)
        pane_copy!(gg[order],i3,1:ws.nstate,g[order],i1,i2,ws.state_index,n)
    end
end

"""
    function make_hh!(hh, g, gg, order, ws)
computes and assembles derivatives of function
    hh(y,u,σ,ϵ) = [y_s; g(y_s,u,σ); g_fwrd(g_state(y_s,u,σ),ϵ,σ); u]
with respect to [y_s, u, σ, ϵ]
"""  
function  make_hh!(hh, g, gg, order, ws)
    if order == 1
        for i = 1:ws.nstate + ws.nshock
            hh[1][i,i] = 1.0
        end
        vh1 = view(hh[1],ws.nstate + (1:ws.nvar),1:(ws.nstate+ws.nshock+1))
        copy!(vh1,g[1])
        n = ws.nstate + 2*ws.nshock + 1
        vh2 = view(hh[1],ws.nstate + ws.nvar + (1:ws.nfwrd),1:n)
        vg2 = view(g[1],ws.fwrd_index,:)
        A_mul_B!(vh2, vg2, gg[1])
        row = ws.nstate + ws.nvar + ws.nfwrd 
        col = ws.nstate
        for i = 1:ws.nshock
            hh[1][row + i, col + i] = 1.0
        end
    else
        # CHECK row order !!!!
        # derivatives of g() for forward looking variables
        copy!(ws.gfwrd[order],view(g[order],ws.fwrd_index,:))
        # derivatives for g(g(y,u,σ),ϵ,σ)
        vh1 = view(hh[order],ws.nstate + ws.ncur + (1:ws.nfwrd),:)
        partial_faa_di_bruno!(vh1, ws.gfwrd, gg, order, ws.faa_di_bruno_ws_1)

        i1 = CartesianIndex(1,(repmat([1], order - 1))...)
        i2 = CartesianIndex(1,(repmat([ws.nstate + ws.nshock + 1], order - 1))...)
        hdims = Tuple(repmat([ws.nstate + 2*ws.nshock + 1],order))
        pane_copy!(hh[order],hdims, ws.cur_index + ws.nstate, g[order], i1, i2, 1:ws.nvar, ws.nstate + ws.nshock + 1)
    end        
end

function pane_copy!(dest,dims,dest_row_range,src,begin_index,end_index,src_row_range,column_nbr)
        r1 = 1:column_nbr
        r2 = r1
        for i in CartesianRange(begin_index,end_index)
            j = sub2ind(dims,(i.I)...) - 1
            v1 = view(dest, dest_row_range, j + r1)
            v2 = view(src, src_row_range, r2)
            copy!(v1,v2)
            r2 += column_nbr
        end
end    

function make_d1!(ws)
    inc1 = ws.nstate
    inc2 = ws.nstate+2*ws.nshock+1
    for j=1:ws.nstate
        col1 = j
        col2 = j
        for k=1:ws.nstate
            for i=1:ws.nvar
                ws.rhs1[i,col1] = -ws.rhs[i,col2]
            end
            col1 += inc1
            col2 += inc2
        end
    end
end

"""
function make_a1!(a::Matrix{Float64}, f::Vector{Matrix{Float64}},
                  g::Vector{Matrix{Float64}}, ncur::Int64,
                  cur_index::Vector{Int64}, nvar::Int64,
                  nstate::Int64, nfwrd::Int64,
                  fwrd_index::Vector{Int64},
                  state_index::Vector{Int64})

updates matrix a with f_0 + f_+g_1 
"""    
function make_a1!(a::Matrix{Float64}, f::Vector{Matrix{Float64}},
                  g::Vector{Matrix{Float64}}, ncur::Int64,
                  cur_index::Vector{Int64}, nvar::Int64,
                  nstate::Int64, nfwrd::Int64,
                  fwrd_index::Vector{Int64},
                  state_index::Vector{Int64})
    
    so = nstate*nvar + 1
    @inbounds for i=1:ncur
        copy!(a,(cur_index[i]-1)*nvar+1,f[1],so,nvar)
        so += nvar
    end
    @inbounds for i = 1:nstate
        for j=1:nstate
            x = 0.0
            @simd for k=1:nfwrd
                x += f[1][j, nstate + ncur + k]*g[1][fwrd_index[k], i]
            end
            a[j,state_index[i]] += x
        end
    end
end

function make_rhs_1_1!(rhs1::AbstractMatrix, rhs::AbstractMatrix, rs::Range{Int64}, rd::Range{Int64}, n::Int64, inc::Int64, order::Int64)
    @inbounds if order > 1
        rs_ = rs
        rd_ = rd
        inc1 = inc^(order-1)
        n1 = n^(order-1)
        for i=1:n
            make_rhs_1_1!(rhs1, rhs, rs_, rd_, n, inc, order - 1)
            rs_ += inc1
            rd_ += n1
        end
    else
        v1 = view(rhs,:,rs)
        v2 = view(rhs1,:,rd)
        v2 .= -v1
    end
end

function make_rhs_1!(rhs1::Matrix{Float64}, rhs::Matrix{Float64}, nstate::Int64,
                     nshock::Int64, nvar::Int64, order::Int64)
    rs = 1:nstate
    rd = 1:nstate
    inc = (nstate + 2*nshock + 1)
    make_rhs_1_1!(rhs1, rhs, rs, rd, nstate, inc, order) 
end

function store_results_1_1!(rhs1::AbstractMatrix, rhs::AbstractMatrix, rs::Range{Int64}, rd::Range{Int64}, n::Int64, inc::Int64, order::Int64)
    @inbounds if order > 1
        rs_ = rs
        rd_ = rd
        inc1 = inc^(order-1)
        n1 = n^(order-1)
        for i = 1:n
            store_results_1_1!(rhs1, rhs, rs_, rd_, n, inc, order - 1)
            rs_ += n1
            rd_ += inc1
        end
    else
        v1 = view(rhs,:,rs)
        v2 = view(rhs1,:,rd)
        v2 .= v1
    end
end

function store_results_1!(result::Matrix{Float64}, gs_ws_result::Matrix{Float64}, nstate::Int64, nshock::Int64, nvar::Int64, order::Int64)
    rs = 1:nstate
    rd = 1:nstate
    inc = (nstate + nshock + 1)
    store_results_1_1!(result, gs_ws_result, rs, rd, nstate, inc, order) 
end

function make_gs_su!(gs_su::Matrix{Float64}, g::Matrix{Float64}, nstate::Int64, nshock::Int64, state_index::Vector{Int64})
    @inbounds for i = 1:(nstate + nshock)
        @simd for j = 1:nstate
            gs_su[j,i] = g[state_index[j],i]
        end
    end
end

function make_gykf_1!(gykf::Matrix{Float64}, g::Matrix{Float64}, rs::Range{Int64}, rd::Range{Int64}, n::Int64, inc::Int64, fwrd_index::Vector{Int64}, order::Int64)
    if order > 1
        rs_ = rs
        rd_ = rd
        inc1 = inc^(order-1)
        n1 = n^(order-1)
        for i = 1:n
            make_gykf_1!(gykf, g, rs, rd, n, inc, fwrd_index, order - 1)
            rs_ += inc1
            rd_ += n1
        end
    else
        v1 = view(g,fwrd_index,rs_)
        v2 = view(gykf,:,rd_)
        v2 .= v1
    end
end

function make_gykf!(gykf::Matrix{Float64}, g::Matrix{Float64}, nstate::Int64, nfwrd::Int64, nshock::Int64, fwrd_index::Vector{Int64}, order::Int64)
    rs = 1:nstate
    rd = 1:nstate
    inc = (nstate + nshock + 1)
    make_gykf_1!(gykf, g, rs, rd, n, inc, fwrd_index, order)
end

function make_rhs_2_1!(rhs1::AbstractMatrix, rhs::AbstractMatrix, os::Int64, od::Int64, n::Int64, inc::Int64, nvar::Int64, blocksize::Int64, order::Int64)
    @inbounds if order > 1
        os_ = os
        od_ = od
        inc1 = nvar*inc^(order-1)
        n1 = nvar*n^(order-1)
        for i= 1:n
            make_rhs_2_1!(rhs1, rhs, os_, od_, n, inc, nvar, blocksize, order - 1)
            os_ += inc1
            od_ += n1
        end
    else
        println(os/nvar + 1)
        k1 = os + 1
        k2 = od + 1
        @simd for i = 1:blocksize
            rhs1[k2] = -rhs1[k2] - rhs[k1]
            k1 += 1
            k2 += 1
        end
    end
end

function make_rhs_2!(rhs1::Matrix{Float64}, rhs::Matrix{Float64}, nstate::Int64,
                     nshock::Int64, nvar::Int64, order::Int64)
    inc = nstate + 2*nshock + 1
    os = nvar*inc^(order-1)
    od = 0
    make_rhs_2_1!(rhs1, rhs, os, od, nstate + nshock, inc, nvar, nvar*(nstate + nshock), order) 
end

function make_rhs_2!(rhs1::Matrix{Float64}, rhs::Matrix{Float64}, nstate::Int64, nshock::Int64, nvar::Int64)
    dcol = 1
    inc = nstate + 2*nshock + 1
    base = nstate*inc + 1
    @inbounds for i=1:nshock
        scol = base 
        for j = 1:(nstate + nshock)
            @simd for k = 1:nvar
                rhs1[k,dcol] = -rhs1[k,dcol] - rhs[k,scol]
            end
            dcol += 1
            scol +=  1
        end
        base += inc
    end
end

"""
    function compute_derivatives_wr_shocks!(ws::KOrderWs,f,g,order)
computes g_su and g_uu
It solves
    (f_+*g_y + f_0)X = -(D + f_+*g_yy*(gu ⊗ [gs gu]) 
"""
function compute_derivatives_wr_shocks!(ws::KOrderWs,f::Vector{Matrix{Float64}},g::Vector{Matrix{Float64}},order::Int64)
    fp = view(f[1],:,ws.nstate + ws.ncur + (1:ws.nfwrd))

    make_gs_su!(ws.gs_su, g[1], ws.nstate, ws.nshock, ws.state_index)
    make_gykf!(ws.gykf, g[order], ws.nstate, ws.nfwrd, ws.nshock, ws.fwrd_index)

    gu = view(ws.gs_su,:,ws.nstate + (1:ws.nshock))
    vrhs1 = view(ws.rhs1,:,1:(ws.nshock*(ws.nstate+ws.nshock)))
    a_mul_b_kron_c_d!(vrhs1,fp,ws.gykf,gu,ws.gs_su,order,ws.work1,ws.work2)

    make_rhs_2!(ws.rhs1, ws.rhs, ws.nstate, ws.nshock, ws.nvar)
    linsolve_core!(ws.linsolve_ws_1,Ref{UInt8}('N'),ws.a,vrhs1)
end
    
function store_results_2!(result::Matrix{Float64}, nstate::Int64, nshock::Int64, nvar::Int64, rhs1::Matrix{Float64})
    soffset = 1
    base1 = nstate*(nstate + nshock + 1)*nvar + 1
    base2 = nstate*nvar + 1
    inc = (nstate + nshock + 1)*nvar
    @inbounds for i=1:nshock
        doffset1 = base1
        doffset2 = base2
        for j=1:(nstate + nshock)
            copy!(result, doffset1, rhs1, soffset, nvar)
            if j <= nstate
                copy!(result, doffset2, rhs1, soffset, nvar)
                doffset2 += inc
            end
            doffset1 += nvar
            soffset +=  nvar
        end
        base1 += (nstate + nshock + 1)*nvar
        base2 += nvar
    end
end

function make_gsk!(g::Vector{Matrix{Float64}},
                   f::Vector{Matrix{Float64}},
                   moments::Vector{Float64}, a::Matrix{Float64},
                   rhs::Matrix{Float64}, rhs1::Matrix{Float64},
                   nfwrd::Int64, nstate::Int64, nvar::Int64,
                   ncur::Int64, nshock::Int64,
                   fwrd_index::Vector{Int64},
                   linsolve_ws_1::LinSolveWS, work1::Vector{Float64},
                   work2::Vector{Float64})
    
    @inbounds for i=1:nfwrd
        @simd for j=1:nvar
            a[j,fwrd_index[i]] += f[1][j, nstate + ncur + i]
        end
    end

    nshock2 = nshock*nshock
    vg = view(work2,1:(nfwrd*nshock2))
    offset = nstate*(nstate+nshock+1) + nstate + 1
    drow = 1
    @inbounds for i=1:nshock
        scol = offset
        for j=1:nshock
            @simd for k = 1:nfwrd
                work2[drow] = g[2][fwrd_index[k],scol]
                drow += 1
            end
            scol += 1
        end
        offset += nstate + nshock + 1
    end
    vfplus = view(f[1],:,nstate + ncur + (1:nfwrd))
    vg1 = reshape(vg,nfwrd,nshock2)
    vwork1 = reshape(view(work1,1:(nvar*nshock2)),nvar,nshock2)
    A_mul_B!(vwork1,vfplus,vg1)
    
    vrhs1 = view(rhs1,:,1:nshock2)
    offset = (nstate + nshock + 1)*(nstate + 2*nshock + 1) + nstate + nshock + 2
    dcol = 1
    @inbounds for i=1:nshock
        scol = offset
        for j=1:nshock
            @simd for k=1:nvar
                vrhs1[k,dcol] = -rhs[k,scol] - vwork1[k,dcol]
            end
            dcol += 1
            scol += 1
        end
        offset += nstate + 2*nshock + 1
    end
    
    vwork2 = view(work2,1:nvar)
    A_mul_B!(vwork2,vrhs1,moments)
    linsolve_core!(linsolve_ws_1,Ref{UInt8}('N'),a,vwork2)
    dcol = (nstate + nshock + 1)
    dcol2 = ((dcol - 1)*dcol + nstate + nshock)*nvar + 1
    copy!(g[2],dcol2,vwork2,1,nvar)
end

"""
    function k_order_solution!(g,f,moments,order,ws)
solves (f^1_0 + f^1_+ gx)X + f^1_+ X (gx ⊗ ... ⊗ gx) = D


"""
function k_order_solution!(g,f,moments,order,ws)
    gg = ws.gg::Vector{Matrix{Float64}}
    hh = ws.hh::Vector{Matrix{Float64}}
    rhs = ws.rhs::Matrix{Float64}
    rhs1 = ws.rhs1::Matrix{Float64}
    faa_di_bruno_ws = ws.faa_di_bruno_ws_2::FaaDiBrunoWs
    nfwrd = ws.nfwrd::Int64
    fwrd_index = ws.fwrd_index::Vector{Int64}
    nstate = ws.nstate::Int64
    state_index = ws.state_index::Vector{Int64}
    ncur = ws.ncur::Int64
    cur_index = ws.cur_index::Vector{Int64}
    nvar = ws.nvar::Int64
    nshock = ws.nshock::Int64
    a = ws.a::Matrix{Float64}
    b = ws.b::Matrix{Float64}
    linsolve_ws = ws.linsolve_ws_1::LinSolveWS
    work1 = ws.work1::Vector{Float64}
    work2 = ws.work2::Vector{Float64}
    gs_ws = ws.gs_ws::EyePlusAtKronBWS
    gs_ws_result = gs_ws.result::Matrix{Float64}
    
    make_gg!(gg, g, order-1, ws)
    make_hh!(hh, g, gg, order-1, ws)
    partial_faa_di_bruno!(rhs,f,hh,order,faa_di_bruno_ws)
    # select only endogenous state variables on the RHS
    # make_d1!(ws)
    if order == 2
        target = copy(a)
        make_a1!(a, f, g, ncur, cur_index, nvar, nstate, nfwrd, fwrd_index, state_index)
        
    end
    @inbounds for i = 1:nfwrd
        col1 = fwrd_index[i]
        col2 = nstate + ncur + i
        @simd for j=1:nvar
            b[j, col1] = f[1][j, col2]
        end
    end
    c = g[1][state_index,1:nstate]
    make_rhs_1!(rhs1, rhs, nstate, nshock, nvar, order)

    d = Vector{Float64}(nvar*nstate^order)
    copy!(d, 1, rhs1, 1, nvar*nstate^order)
    println("generalized_sylvester")
    generalized_sylvester_solver!(a,b,c,d,order,gs_ws)
    println("end generalized_sylvester")
    store_results_1!(g[order], gs_ws_result, nstate, nshock, nvar, order)
    compute_derivatives_wr_shocks!(ws,f,g,order)
    store_results_2!(g[order], nstate, nshock, nvar, rhs1, order)
    make_gsk!(g, f, moments[2], a, rhs, rhs1,
              nfwrd, nstate, nvar, ncur, nshock,
              fwrd_index, linsolve_ws, work1, work2)
end

end
