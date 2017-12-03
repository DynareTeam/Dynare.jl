module KOrder

using Base.Test
using ...DynLinAlg.Kronecker
using ...DynLinAlg.linsolve_algo
import Base.LinAlg.BLAS: gemm!
import ..GeneralizedSylvester: EyePlusAtKronBWS, generalized_sylvester_solver!
import ...FaaDiBruno: partial_faa_di_bruno!, FaaDiBrunoWs
export make_gg!, make_hh!, k_order_solution!, KOrderWs

struct KOrderWs
    nvar::Integer
    nfwrd::Integer
    nstate::Integer
    ncur::Integer
    nshock::Integer
    fwrd_index::Array{Integer}
    state_index::Array{Integer}
    cur_index::Array{Integer}
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
        new(nvar,nfwrd,nstate,ncur,nshock,fwrd_index,state_index,cur_index,state_range,gfwrd,gg,hh,
            rhs,rhs1,gykf,gs_su,a,b,work1,work2,faa_di_bruno_ws_1,faa_di_bruno_ws_2,linsolve_ws_1)
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
    hh(y,u,σ,ϵ) = [g_fwrd(g_state(y_s,u,σ),ϵ,σ);g(y_s,u,σ);y_s;u]
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
    else
        # CHECK row order !!!!
        # derivatives of g() for forward looking variables
        copy!(ws.gfwrd[order],view(g[order],ws.fwrd_index,:))
        # derivatives for g(g(y,u,σ),ϵ,σ)
        vh1 = view(hh[order],1:ws.nfwrd,:)
        # partial_faa_di_bruno needs Matrix argument! TO BE FIXED
        tmp = copy(vh1)
        partial_faa_di_bruno!(tmp,ws.gfwrd,gg,order,ws.faa_di_bruno_ws_1)
        copy!(vh1,tmp)

        i1 = CartesianIndex(1,(repmat([1], order - 1))...)
        i2 = CartesianIndex(1,(repmat([ws.nstate + ws.nshock + 1], order - 1))...)
        hdims = Tuple(repmat([ws.nstate + 2*ws.nshock + 1],order))
        pane_copy!(hh[order],hdims, ws.cur_index + ws.nfwrd, g[order], i1, i2, 1:ws.nvar, ws.nstate + ws.nshock + 1)
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

function make_a1!(ws,f,g)
    so = ws.nfwrd*ws.nvar + 1
    for i=1:ws.ncur
        copy!(ws.a,(ws.cur_index[i]-1)*ws.nvar+1,f[1],so,ws.nvar)
        so += ws.nvar
    end
    for i = 1:ws.nstate
        for j=1:ws.nstate
            x = 0.0
            for k=1:ws.nfwrd
                x += f[1][j,ws.nstate + ws.ncur + k]*g[1][ws.fwrd_index[k],i]
            end
            ws.a[j,ws.state_index[i]] += x
        end
    end
end

function make_rhs_1(ws::KOrderWs)
    dcol = 1
    base = 1
    inc = (ws.nstate + 2*ws.nshock + 1)
    for i=1:ws.nstate
        scol = base 
        for j=1:ws.nstate
            vrhs1 = view(ws.rhs1,:,dcol)
            vrhs = view(ws.rhs,:,scol)
            # to be checked for optimization
            vrhs1 = -vrhs
            dcol += 1
            scol += 1
        end
        base += inc
    end
    view(ws.rhs1,:,1:(ws.nstate*ws.nstate))
end

function store_results_1!(result,ws)
    soffset = 1
    base = 1
    inc = (ws.nstate + ws.nshock + 1)*ws.nvar
    for i=1:ws.nstate
        doffset = base 
        for j=1:ws.nstate
            copy!(result, doffset, ws.rhs1, soffset, ws.nvar)
            doffset += ws.nvar
            soffset +=  ws.nvar
        end
        base += inc
    end
end

"""
    function make_rhs2!(ws::KOrderWs,f,g,order)
computes g_su and g_uu
It solves
    (f_+*g_y + f_0)X = -(D + f_+*g_yy*(gu ⊗ [gs gu]) 
"""
function make_rhs2!(ws::KOrderWs,f,g,order)
    fp = view(f[1],:,ws.nstate + ws.ncur + (1:ws.nfwrd))

    for i = 1:(ws.nstate + ws.nshock)
        for j = 1:ws.nstate
            ws.gs_su[j,i] = g[1][ws.state_index[j],i]
        end
    end
    
    dcol = 1
    base = 1
    for i = 1:ws.nstate
        scol = base
        for j = 1:ws.nstate
            for k=1:ws.nfwrd
                ws.gykf[k,dcol] = g[order][ws.fwrd_index[k],scol]
            end
            dcol += 1
            scol += 1
        end
        base += ws.nstate + ws.nshock + 1
    end

    gu = view(ws.gs_su,:,ws.nstate + (1:ws.nshock))
    vrhs1 = view(ws.rhs1,:,1:(ws.nshock*(ws.nstate+ws.nshock)))
    a_mul_b_kron_c_d!(vrhs1,fp,ws.gykf,gu,ws.gs_su,order,ws.work1,ws.work2)
    
    dcol = 1
    inc = ws.nstate + 2*ws.nshock + 1
    base = ws.nstate*inc + 1
    for i=1:ws.nshock
        scol = base 
        for j = 1:(ws.nstate + ws.nshock)
            for k = 1:ws.nvar
                ws.rhs1[k,dcol] = -ws.rhs1[k,dcol] - ws.rhs[k,scol]
            end
            dcol += 1
            scol +=  1
        end
        base += inc
    end
    linsolve_core!(ws.linsolve_ws_1,Ref{UInt8}('N'),ws.a,vrhs1)
end
    
function store_results_2!(result,ws)
    soffset = 1
    base1 = ws.nstate*(ws.nstate + ws.nshock + 1)*ws.nvar + 1
    base2 = ws.nstate*ws.nvar + 1
    inc = (ws.nstate + ws.nshock + 1)*ws.nvar
    for i=1:ws.nshock
        doffset1 = base1
        doffset2 = base2
        for j=1:(ws.nstate + ws.nshock)
            copy!(result, doffset1, ws.rhs1, soffset, ws.nvar)
            if j <= ws.nstate
                copy!(result, doffset2, ws.rhs1, soffset, ws.nvar)
                doffset2 += inc
            end
            doffset1 += ws.nvar
            soffset +=  ws.nvar
        end
        base1 += (ws.nstate + ws.nshock + 1)*ws.nvar
        base2 += ws.nvar
    end
end

function make_gsk!(ws,f,g,moments)
    for i=1:ws.nfwrd
        for j=1:ws.nvar
            ws.a[j,ws.fwrd_index[i]] += f[1][j, ws.nstate + ws.ncur + i]
        end
    end

    nshock2 = ws.nshock*ws.nshock
    vg = view(ws.work2,1:(ws.nfwrd*nshock2))
    offset = ws.nstate*(ws.nstate+ws.nshock+1) + ws.nstate + 1
    drow = 1
    for i=1:ws.nshock
        scol = offset
        for j=1:ws.nshock
            for k = 1:ws.nfwrd
                vg[drow] = g[2][ws.fwrd_index[k],scol]
                drow += 1
            end
            scol += 1
        end
        offset += ws.nstate + ws.nshock + 1
    end
    vfplus = view(f[1],:,ws.nstate + ws.ncur + (1:ws.nfwrd))
    vg1 = reshape(vg,ws.nfwrd,nshock2)
    vwork1 = reshape(view(ws.work1,1:(ws.nvar*nshock2)),ws.nvar,nshock2)
    A_mul_B!(vwork1,vfplus,vg1)
    
    vrhs1 = view(ws.rhs1,:,1:nshock2)
    offset = (ws.nstate + ws.nshock + 1)*(ws.nstate + 2*ws.nshock + 1) + ws.nstate + ws.nshock + 2
    dcol = 1
    for i=1:ws.nshock
        scol = offset
        for j=1:ws.nshock
            for k=1:ws.nvar
                vrhs1[k,dcol] = -ws.rhs[k,scol] - vwork1[k,dcol]
            end
            dcol += 1
            scol += 1
        end
        offset += ws.nstate + 2*ws.nshock + 1
    end
    
    vwork2 = view(ws.work2,1:ws.nvar)
    A_mul_B!(vwork2,vrhs1,moments)
    linsolve_core!(ws.linsolve_ws_1,Ref{UInt8}('N'),ws.a,vwork2)
    dcol = (ws.nstate + ws.nshock + 1)
    dcol2 = ((dcol - 1)*dcol + ws.nstate + ws.nshock)*ws.nvar + 1
    copy!(g[2],dcol2,vwork2,1,ws.nvar)
end

"""
    function k_order_solution!(g,f,moments,order,ws)
solves (f^1_0 + f^1_+ gx)X + f^1_+ X (gx ⊗ ... ⊗ gx) = D


"""
function k_order_solution!(g,f,moments,order,ws)
    @time make_gg!(ws.gg, g, order-1, ws)
    @time make_hh!(ws.hh, g, ws.gg, order-1, ws)
    @time partial_faa_di_bruno!(ws.rhs,f,ws.hh,order,ws.faa_di_bruno_ws_2)
    # select only endogenous state variables on the RHS
    make_d1!(ws)
    if order == 2
        @time make_a1!(ws,f,g)
    end
    for i = 1:ws.nfwrd
        col1 = ws.fwrd_index[i]
        col2 = ws.nstate + ws.ncur + i
        for j=1:ws.nvar
            ws.b[j, col1] = f[1][j, col2]
        end
    end
    c = g[1][ws.state_index,1:ws.nstate]
    @time gs_ws = EyePlusAtKronBWS(ws.nvar,ws.nvar,ws.nstate,order)
    rhs1 = make_rhs_1(ws)

    @time generalized_sylvester_solver!(ws.a,ws.b,c,rhs1,order,gs_ws)
    @time store_results_1!(g[2],ws)
    @time make_rhs2!(ws,f,g,order)
    @time store_results_2!(g[2],ws)

    @time make_gsk!(ws,f,g,moments[2])
end

end
