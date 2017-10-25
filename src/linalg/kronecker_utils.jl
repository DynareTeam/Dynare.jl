module Kronecker

using Base.Test
import Base.convert
export a_mul_kron_b!, a_mul_b_kron_c!, kron_at_kron_b_mul_c!, a_mul_b_kron_c_d!

"""
    a_mul_kron_b!(c::AbstractVector, a::AbstractVecOrMat, b::AbstractMatrix, order::Int64)

Performs a*(b ⊗ b ⊗ ... ⊗ b). The solution is returned in matrix c. order indicates the number of occurences of b

We use vec(a*(b ⊗ b ⊗ ... ⊗ b)) = (b' ⊗ b' ⊗ ... ⊗ b' \otimes I)vec(a)

"""
function a_mul_kron_b!(c::AbstractMatrix, a::AbstractMatrix, b::AbstractMatrix, order::Int64, work1::AbstractVector, work2::AbstractVector)
    ma, na = size(a)
    mb, nb = size(b)
    mc, nc = size(c)
    mborder = mb^order
    nborder = nb^order
    na == mborder || throw(DimensionMismatch("The number of columns of a, $na, doesn't match the number of rows of b, $mb, times order = $order"))
    mc == ma || throw(DimensionMismatch("The number of rows of c, $mc, doesn't match the number of rows of a, $ma"))
    nc == nborder || throw(DimensionMismatch("The number of columns of c, $nc, doesn't match the number of columns of b, $nb, times order = $order"))

    # copy input
    copy!(work1,a)
    n1 = ma*na
    n2 = Int(n1*nb/mb)
    v1 = view(work1,1:n1)
    v2 = view(work2,1:n2)
    for q=0:order-1
        kron_mul_elem_t!(v2,b,v1,mb^(order-q-1),nb^q*ma)
        if q < order - 1
            # update dimension
            n2 = Int(n2*nb/mb)
            # swap and resize work vectors
            v1parent = v1.parent
            v1 = v2
            v2 = view(v1parent,1:n2)
        end
    end
    copy!(c,v2)
end
    
function a_mul_kron_b!(c::AbstractMatrix, a::AbstractMatrix, b::AbstractMatrix, order::Int64)
    ma, na = size(a)
    mb, nb = size(b)
    mc, nc = size(c)
    mborder = mb^order
    nborder = nb^order
    na == mborder || throw(DimensionMismatch("The number of columns of a, $na, doesn't match the number of rows of b, $mb, times order = $order"))
    mc == ma || throw(DimensionMismatch("The number of rows of c, $mc, doesn't match the number of rows of a, $ma"))
    nc == nborder || throw(DimensionMismatch("The number of columns of c, $nc, doesn't match the number of columns of b, $nb, times order = $order"))
    mb == nb || throw(DimensionMismatch("B must be a square matrix"))
 
    avec = vec(a)
    cvec = vec(c)
    for q=0:order-1
        vavec = view(avec,1:ma*mb^(order-q)*nb^q)
        vcvec = view(cvec,1:ma*mb^(order-q-1)*nb^(q+1))
        kron_mul_elem_t!(vcvec,b,vavec,mb^(order-q-1),nb^q*ma)
        if q < order - 1
            copy!(avec,vcvec)
        end
    end
end
    
"""
    a_mul_b_kron_c!(d::AbstractVecOrMat, a::AbstractVecOrMat, b::AbstractMatrix, c::AbstractMatrix, order::Int64)

Performs a*B*(c ⊗ c ⊗ ... ⊗ c). The solution is returned in matrix or vector d. order indicates the number of occurences of c. c must be a square matrix

We use vec(a*b*(c ⊗ c ⊗ ... ⊗ c)) = (c' ⊗ c' ⊗ ... ⊗ c' ⊗ a)vec(b)

"""
function a_mul_b_kron_c!(d::AbstractMatrix, a::AbstractMatrix, b::AbstractMatrix, c::AbstractMatrix, order::Int64)
    ma, na = size(a)
    mb, nb = size(b)
    mc, nc = size(c)
    md, nd = size(d)
    ma == na || throw(DimensionMismatch("a must be a square matrix"))
    mc == nc || throw(DimensionMismatch("c must be a square matrix"))
    na == mb || throw(DimensionMismatch("The number of columns of a, $(size(a,2)), doesn't match the number of rows of b, $(size(b,1))"))
    nb == mc^order || throw(DimensionMismatch("The number of columns of b, $(size(b,2)), doesn't match the number of rows of c, $(size(c,1)), times order, $order"))
    (ma == md && nc^order == nd) || throw(DimensionMismatch("Dimension mismatch for D: $(size(d)) while ($ma, $(nc^order)) was expected"))
    A_mul_B!(d,a,b)
    copy!(b,d)
    bvec = vec(b)
    dvec = vec(d)
    for q=0:order-1
        kron_mul_elem_t!(dvec,c,bvec,mc^(order-q-1),nc^q*ma)
        if q < order - 1
            copy!(bvec,dvec)
        end
    end
end

"""
    function kron_at_kron_b_mul_c!(a::AbstractMatrix, order::Int64, b::AbstractMatrix, c::AbstractVector, work::AbstractVector)
updates c by (a^T ⊗ a^T ⊗ ... ⊗ a^T ⊗ b)c
""" 
function kron_at_kron_b_mul_c!(a::AbstractMatrix, order::Int64, b::AbstractMatrix, c::AbstractVector, work::AbstractVector)
    copy!(work,c)
    ma,na = size(a)
    mb,nb = size(b)
    maorder = ma^order
    naorder = na^order
    length(work) == maorder*nb  || throw(DimensionMismatch("The dimension of vector , $(length(b)) doesn't correspond to order, ($p, $q)  and the dimension of the matrices a, $(size(a)), and b, $(size(b))"))
    work1 = convert(Matrix{Float64},reshape(work,nb,naorder))
    c1 = convert(Matrix{Float64},reshape(c,nb,maorder))
    A_mul_B!(work1,b,c1)
    copy!(c1,work1)
    for q = 0:order-1
        kron_mul_elem_t!(work,a,c,mb^(order-q-1),nb^q*nb)
        copy!(c,work)
    end
end

function a_mul_kron_b!(c::AbstractMatrix, a::AbstractMatrix, b::AbstractVector,work::AbstractVector)
    ma, na = size(a)
    order = length(b)
    mc, nc = size(c)
    mborder = 1
    nborder = 1
    for i = 1:order
        mb, nb = size(b[i])
        mborder *= mb
        nborder *= nb
    end
    na == mborder || throw(DimensionMismatch("The number of columns of a, $na, doesn't match the number of rows of matrices in b, $mborder"))
    mc == ma || throw(DimensionMismatch("The number of rows of c, $mc, doesn't match the number of rows of a, $ma"))
    nc == nborder || throw(DimensionMismatch("The number of columns of c, $nc, doesn't match the number of columns of matrices in b, $nborder"))
    mborder <= nborder || throw(DimensionMismatch("the product of the number of rows of the b matrices needs to be smaller or equal to the product of the number of columns. Otherwise, you need to call a_mul_kron_b!(c::AbstractMatrix, a::AbstractMatrix, b::Vector{AbstractMatrix},work::AbstractVector)"))
 
    copy!(work,a)
    mb, nb = size(b[1])
    vwork = view(work,1:ma*na)
    cvec = vec(c)
    mc = ma*na
    p = Int(mborder/size(b[order],1))
    q = ma
    for i = order:-1:1
        mb, nb = size(b[i])
        mc = Int(mc*nb/mb)
        vcvec = view(cvec,1:mc)
        kron_mul_elem_t!(vcvec,b[i],vwork,p,q)
        if i > 1
            p = Int(p/mb)
            q = q*nb
            vwork = view(work,1:mc)
            copy!(vwork,vcvec)
        end
    end
end
    
"""
    a_mul_b_kron_c_d!(d::AbstractVecOrMat, a::AbstractVecOrMat, b::AbstractMatrix, c::AbstractMatrix, order::Int64)

Performs a*b*(c ⊗ d ⊗ ... ⊗ d). The solution is returned in matrix or vector e order indicates the number of occurences of d. c and d must be square matrices

We use vec(a*b*(c ⊗ d ⊗ ... ⊗ d)) = (c' ⊗ d' ⊗ ... ⊗ d' ⊗ a)vec(b)

"""
function a_mul_b_kron_c_d!(e::AbstractMatrix, a::AbstractMatrix, b::AbstractMatrix, c::AbstractMatrix, d::AbstractMatrix, order::Int64, work1::AbstractVector, work2::AbstractVector)

    ma, na = size(a)
    mb, nb = size(b)
    mc, nc = size(c)
    md, nd = size(d)
    me, ne = size(e)
    na == mb || throw(DimensionMismatch("The number of columns of a, $(size(a,2)), doesn't match the number of rows of b, $(size(b,1))"))
    nb == mc*md^order || throw(DimensionMismatch("The number of columns of b, $(size(b,2)), doesn't match the number of rows of c, $(size(c,1)), times order, $order"))
    (ma == me && nc*nd^order == ne) || throw(DimensionMismatch("Dimension mismatch for e: $(size(e)) while ($ma, $(nc*nd^order)) was expected"))
    v1 = reshape(view(work1,1:ma*nb),ma,nb)
    A_mul_B!(v1,a,b)
    copy!(work2,work1)
    vw2 = vec(work2)
    for q=0:order-1
        v2 = view(vw2,1:ma*mc*nd^(order-q)*nd^q)
        v1 = view(work1,1:ma*mc*md^(order-q)*nd^q)
        kron_mul_elem_t!(v2,d,v1,mc*md^(order-q-1),nd^q*ma)
        copy!(work1,v2)
    end
    v2 = view(work2,1:ma*nc*nd^order)
    v1 = view(work1,1:ma*mc*nd^order)
    kron_mul_elem_t!(v2,c,v1,1,nd^order*ma)
    copy!(e,v2)
end

convert(::Type{Array{Float64, 2}}, x::Base.ReshapedArray{Float64,2,SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},Tuple{}}) = unsafe_wrap(Array,pointer(x.parent.parent,x.parent.indexes[1][1]),x.dims)


"""
    kron_mul_elem!(c::AbstractVector, a::AbstractMatrix, b::AbstractVector, p::Int64, q::Int64, m::Int64)

Performs (I_n^p \otimes a \otimes I_{m*n^q}) b, where n = size(a,2) and m, an arbitrary constant. The result is stored in c.
"""
function kron_mul_elem!(c::AbstractVector, a::AbstractMatrix, b::AbstractVector, p::Int64, q::Int64, m::Int64)
    ma, na = size(a)
    length(b) == m*ma^p*na^(q+1) || throw(DimensionMismatch("The dimension of vector b, $(length(b)) doesn't correspond to order, ($p, $q)  and the dimensions of the matrix, $(size(a))"))
    length(c) == m*ma^(p+1)*na^q || throw(DimensionMismatch("The dimension of the vector c, $(length(c)) doesn't correspond to order, ($p, $q)  and the dimensions of the matrix, $(size(a))"))


    @inbounds begin
        if m == 1 && p + q == 0
            # a*b
            A_mul_B!(c,a,b)
        elseif m == 1 && q == 0
            #  (I_n^p ⊗ a)*b = vec(a*[b_1 b_2 ... b_p])
            b = convert(Array{Float64,2},reshape(b,na,ma^p))
            c = convert(Array{Float64,2},reshape(c,ma,ma^p))
            A_mul_B!(c,a,b)
        elseif p == 0
            # (a ⊗ I_{m*n^q})*b = (b'*(a' ⊗ I_{m*n^q}))' = vec(reshape(b,m*n^q,n)*a')
            @time b = convert(Array{Float64,2},reshape(b,m*ma^p*na^q,na))
            @time c = convert(Array{Float64,2},reshape(c,m*ma^p*na^q,ma))
            @time A_mul_Bt!(c,b,a)
        else
            # (I_{n^p} ⊗ a ⊗ I_{m*n^q})*b = vec([(a ⊗ I_{m*n^q})*b_1 (a ⊗ I_{m*n^q})*b_2 ... (a ⊗ I_{m*n^q})*b_{n^p}])
            mnq = m*na^q
            mnq1 = mnq*na
            qrangeb = 1:mnq1
            mnq2 = mnq*ma
            qrangec = 1:mnq2
            for i=1:ma^p
                bi = convert(Array{Float64,2},reshape(view(b,qrangeb),m*na^q,na))
                ci = convert(Array{Float64,2},reshape(view(c,qrangec),m*na^q,ma))
                # (a ⊗ I_{m*n^q})*b = (b'*(a' ⊗ I_{m*n^q}))' = vec(reshape(b',m*n^q,n)*a')
                A_mul_Bt!(ci,bi,a)
                qrangeb += mnq1
                qrangec += mnq2
            end
        end
    end
end

"""
    kron_mul_elem_t!(p::Int64, q::Int64, m::Int64, a::AbstractMatrix, b::AbstractVector, c::AbstractVector)

Performs (I_p \otimes a' \otimes I_q) b, where m,n = size(a). The result is stored in c.
"""
function kron_mul_elem_t!(c::AbstractVector, a::AbstractMatrix, b::AbstractVector, p::Int64, q::Int64)
    m, n = size(a)
    length(b) == m*p*q || throw(DimensionMismatch("The dimension of vector b, $(length(b)) doesn't correspond to order, ($p, $q)  and the dimensions of the matrix, $(size(a))"))
    length(c) == n*p*q || throw(DimensionMismatch("The dimension of the vector c, $(length(c)) doesn't correspond to order, ($p, $q)  and the dimensions of the matrix, $(size(a))"))

    @inbounds begin
        if p == 1 && q == 1
            # a'*b
            At_mul_B!(c,a,b)
        elseif q == 1
            #  (I_p ⊗ a')*b = vec(a'*[b_1 b_2 ... b_p])
            b = convert(Array{Float64,2},reshape(b,m,p))
            c = convert(Array{Float64,2},reshape(c,n,p))
            At_mul_B!(c,a,b)
        elseif p == 1
            # (a' ⊗ I_q)*b = (b'*(a ⊗ I_q))' = vec(reshape(b,q,m)*a)
            b = convert(Array{Float64,2},reshape(b,q,m))
            c = convert(Array{Float64,2},reshape(c,q,n))
            A_mul_B!(c,b,a)
        else
            # (I_p ⊗ a' ⊗ I_q)*b = vec([(a' ⊗ I_q)*b_1 (a' ⊗ I_q)*b_2 ... (a' ⊗ I_q)*b_p])
            mq = m*q
            nq = n*q
            qrange1 = 1:mq
            qrange2 = 1:nq
            for i=1:p
                bi = convert(Array{Float64,2},reshape(view(b,qrange1),q,m))
                ci = convert(Array{Float64,2},reshape(view(c,qrange2),q,n))
                # (a' ⊗ I_q)*bi = (bi'*(a ⊗ I_q))' = vec(reshape(bi,q,n)*a)
                A_mul_B!(ci,bi,a)
                qrange1 += mq
                qrange2 += nq
            end
        end
    end
end

end
