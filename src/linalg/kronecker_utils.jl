module Kronecker

using Base.Test
import Base.convert
export a_mul_kron_b!, a_mul_b_kron_c!, kron_at_kron_b_mul_c!, a_mul_b_kron_c_d!, at_mul_b_kron_c!, a_mul_b_kron_ct!
import Base.LinAlg: A_mul_B!, At_mul_B!, A_mul_Bt!, At_mul_Bt
import Base.BLAS: gemm!

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
#        vavec = view(avec,1:ma*mb^(order-q)*nb^q)
#        vcvec = view(cvec,1:ma*mb^(order-q-1)*nb^(q+1))
        kron_mul_elem_t!(cvec,b,avec,mb^(order-q-1),nb^q*ma)
        if q < order - 1
            copy!(avec,cvec)
        end
    end
end
    
"""
    a_mul_b_kron_c!(d::AbstractVecOrMat, a::AbstractVecOrMat, b::AbstractMatrix, c::AbstractMatrix, order::Int64)

Performs a*B*(c ⊗ c ⊗ ... ⊗ c). The solution is returned in matrix or vector d. order indicates the number of occurences of c. c must be a square matrix

We use vec(a*b*(c ⊗ c ⊗ ... ⊗ c)) = (c' ⊗ c' ⊗ ... ⊗ c' ⊗ a)vec(b)

"""
function a_mul_b_kron_c!(d::AbstractMatrix, a::AbstractMatrix, b::AbstractMatrix, c::AbstractMatrix, order::Int64,work1::AbstractVector,work2::AbstractVector)
    ma, na = size(a)
    mb, nb = size(b)
    mc, nc = size(c)
    md, nd = size(d)
    na == mb || throw(DimensionMismatch("The number of columns of a, $(size(a,2)), doesn't match the number of rows of b, $(size(b,1))"))
    nb == mc^order || throw(DimensionMismatch("The number of columns of b, $(size(b,2)), doesn't match the number of rows of c, $(size(c,1)), times order, $order"))
    (ma == md && nc^order == nd) || throw(DimensionMismatch("Dimension mismatch for D: $(size(d)) while ($ma, $(nc^order)) was expected"))
    A_mul_B!(reshape(view(work1,1:ma*nb),ma,nb),a,b)
    for q=0:order-1
        kron_mul_elem_t!(work2,c,work1,mc^(order-q-1),nc^q*ma)
        if q < order - 1
            copy!(work1,work2)
        end
    end
    copy!(d,work2)
end

"""
    function kron_at_kron_b_mul_c!(d::AbstractVector, offset_c::Int64, a::AbstractMatrix, order::Int64, b::AbstractMatrix, c::AbstractVector, offset_cc::Int64, work1::AbstractVector, work2::AbstractVector)
computes d = (a^T ⊗ a^T ⊗ ... ⊗ a^T ⊗ b)c using work vectors work1 and work2
""" 
function kron_at_kron_b_mul_c!(d::AbstractVector, offset_d::Int64, a::AbstractMatrix, order::Int64, b::AbstractMatrix, c::AbstractVector, offset_c::Int64, work1::AbstractVector, work2::AbstractVector, offset_w::Int64)
    mb,nb = size(b)
    if order == 0
        A_mul_B!(d,offset_d,b,1,mb,nb,c,offset_c,1)
    else
        ma, na = size(a)
        #    length(work) == naorder*mb  || throw(DimensionMismatch("The dimension of vector , $(length(c)) doesn't correspond to order, ($order)  and the dimension of the matrices a, $(size(a)), and b, $(size(b))"))
        p = ma^order
        kron_mul_elem!(work1, 1, b, c, offset_c, p, 1)
        p = Int(p/ma)
        q = mb
        for i = 1:order
            kron_mul_elem_t!(work2, offset_w, a, work1, 1, p, q)
            if i < order
                copy!(work1, 1, work2, offset_w, p*na*q)
                p = Int(p/ma)
                q *= na
            end
        end
        copy!(d, offset_d, work2, offset_w, p*na*q)
    end
end

"""
    function kron_at_kron_b_mul_c!(d::AbstractVector, a::AbstractMatrix, order::Int64, b::AbstractMatrix, c::AbstractVector, work1::AbstractVector, work2::AbstractVector)
computes d = (a^T ⊗ a^T ⊗ ... ⊗ a^T ⊗ b)c using work vectors work1 and work2
""" 
function kron_at_kron_b_mul_c!(d::AbstractVector, a::AbstractMatrix, order::Int64, b::AbstractMatrix, c::AbstractVector, work1::AbstractVector, work2::AbstractVector)
    mb,nb = size(b)
    if order == 0
        A_mul_B!(d,1,b,1,mb,nb,c,1,1)
    else
        ma, na = size(a)
        #    length(work) == naorder*mb  || throw(DimensionMismatch("The dimension of vector , $(length(c)) doesn't correspond to order, ($order)  and the dimension of the matrices a, $(size(a)), and b, $(size(b))"))
        p = ma^order
        kron_mul_elem!(work1, b, c, p, 1)
        p = Int(p/ma)
        q = mb
        for i = 1:order
            kron_mul_elem_t!(work2, a, work1, p, q)
            if i < order
                copy!(work1, 1, work2, 1, p*na*q)
                p = Int(p/ma)
                q *= na
            end
        end
        copy!(d,1,work2,1,p*na*q)
    end
end

function kron_at_kron_b_mul_c!(a::AbstractMatrix, order::Int64, b::AbstractMatrix, c::AbstractVector, offset_c::Int64, work::AbstractVector)
    kron_at_kron_b_mul_c!(work, 1, a, order, b, c, offset_c, work, c, offset_c)
    copy!(c, work)
end

function kron_at_kron_b_mul_c!(a::AbstractMatrix, order::Int64, b::AbstractMatrix, c::AbstractVector, work::AbstractVector)
    kron_at_kron_b_mul_c!(work, a, order, b, c, work, c)
    copy!(c, work)
end

function kron_at_mul_b!(c::AbstractVector, a::AbstractMatrix, order::Int64, b::AbstractVector, q::Int64, work1::AbstractVector, work2::AbstractVector)
    ma,na = size(a)
#    length(work) == naorder*mb  || throw(DimensionMismatch("The dimension of vector , $(length(c)) doesn't correspond to order, ($order)  and the dimension of the matrices a, $(size(a)), and b, $(size(b))"))
    p = ma^(order-1)
    s = p*q*na
    copy!(work1,b)
    for i = 1:order
        kron_mul_elem_t!(work2, a, work1, p, q)
        if i < order
            copy!(work1, 1, work2, 1, s)
            p = Int(p/ma)
            q *= na
            s = p*q*na
        end
    end
    copy!(c,1,work2,1,s)
end

function kron_a_mul_b!(c::AbstractVector, a::AbstractMatrix, order::Int64, b::AbstractVector, q::Int64, work1::AbstractVector, work2::AbstractVector)
    ma,na = size(a)
#    length(work) == naorder*mb  || throw(DimensionMismatch("The dimension of vector , $(length(c)) doesn't correspond to order, ($order)  and the dimension of the matrices a, $(size(a)), and b, $(size(b))"))
    p = na^(order-1)
    s = p*q*ma
    copy!(work1,b)
    for i = 1:order
        kron_mul_elem!(work2,a,work1,p,q)
        if i < order
            copy!(work1, 1, work2, 1, s)
            p = Int(p/na)
            q *= ma
            s = p*q*ma
        end
    end
    copy!(c,1,work2,1,s)
end

function at_mul_b_kron_c!(d::AbstractMatrix, a::AbstractMatrix, b::AbstractMatrix, c::AbstractMatrix, order::Int64, work1::AbstractVector, work2::AbstractVector)
    ma, na = size(a)
    mb, nb = size(b)
    mc, nc = size(c)
    if mc <= nc
        At_mul_B!(work1, 1, a, 1, ma, na, b, 1, nb)
        kron_at_mul_b!(vec(d), c, order, work1, na, work1, work2)
    else
        kron_at_mul_b!(work1, c, order, b, na, work1, work2)
        At_mul_B!(vec(d), 1, a, 1, ma, na, work1, 1, nc^order)
    end
end
                      
function a_mul_b_kron_ct!(d::AbstractMatrix, a::AbstractMatrix, b::AbstractMatrix, c::AbstractMatrix, order::Int64, work1::AbstractVector, work2::AbstractVector)
    ma, na = size(a)
    mb, nb = size(b)
    mc, nc = size(c)
    if mc <= nc
        A_mul_B!(work1, 1, a, 1, ma, na, b, 1, nb)
        kron_a_mul_b!(vec(d), c, order, work1, ma, work1, work2)
    else
        kron_a_mul_b!(work1, b, order, c, mb, work1, work2)
        At_mul_B!(vec(d), 1, a, 1, ma, na, work1, 1, nc^order)
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
        kron_mul_elem_t!(cvec,b[i],work,p,q)
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
     nb == mc*md^(order-1) || throw(DimensionMismatch("The number of columns of b, $(size(b,2)), doesn't match the number of rows of c, $(size(c,1)), and d, $(size(d,1)) for order, $order"))
     (ma == me && nc*nd^(order-1) == ne) || throw(DimensionMismatch("Dimension mismatch for e: $(size(e)) while ($ma, $(nc*nd^(order-1))) was expected"))
     A_mul_B!(work1, 1, a, 1, ma, mb, vec(b), 1, nb)
     p = mc*md^(order - 2)
     q = ma
     for i = 0:order - 2
        @time kron_mul_elem_t!(work2, 1, d, work1, 1, p, q)
        copy!(work1,work2)
        p = Int(p/md)
        q *= nd
    end
    @time kron_mul_elem_t!(work2, c, work1, 1,q)
    copy!(e, 1, work2, 1, ma*nc*nd^(order-1))
end

convert(::Type{Array{Float64, 2}}, x::Base.ReshapedArray{Float64,2,SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},Tuple{}}) = unsafe_wrap(Array,pointer(x.parent.parent,x.parent.indexes[1][1]),x.dims)


"""
    kron_mul_elem_t!(p::Int64, q::Int64, m::Int64, a::AbstractMatrix, b::AbstractVector, c::AbstractVector)

Performs (I_p \otimes a \otimes I_q) b, where m,n = size(a). The result is stored in c.
"""
function kron_mul_elem!(c::AbstractVector, offset_c::Int64, a::AbstractMatrix, b::AbstractVector, offset_b::Int64, p::Int64, q::Int64)
    m, n = size(a)
    length(b) >= n*p*q || throw(DimensionMismatch("The dimension of vector b, $(length(b)) doesn't correspond to order, ($p, $q)  and the dimensions of the matrix, $(size(a))"))
    length(c) >= m*p*q || throw(DimensionMismatch("The dimension of the vector c, $(length(c)) doesn't correspond to order, ($p, $q)  and the dimensions of the matrix, $(size(a))"))

    begin
        if p == 1 && q == 1
            # a*b
            A_mul_B!(c, offset_c, a, 1, m, n, b, offset_b, 1)
        elseif q == 1
            #  (I_p ⊗ a)*b = vec(a*[b_1 b_2 ... b_p])
            A_mul_B!(c, offset_c, a, 1, m, n, b, offset_b, p)
        elseif p == 1
            # (a ⊗ I_q)*b = (b'*(a' ⊗ I_q))' = vec(reshape(b,q,m)*a')
            A_mul_Bt!(c, offset_c, b, offset_b, q, n, a, 1, m)
        else
            # (I_p ⊗ a ⊗ I_q)*b = vec([(a ⊗ I_q)*b_1 (a ⊗ I_q)*b_2 ... (a ⊗ I_q)*b_p])
            mq = m*q
            nq = n*q
            for i=1:p
                A_mul_Bt!(c, offset_c, b, offset_b, q, n, a, 1, m)
                offset_b += nq
                offset_c += mq
            end
        end
    end
end


"""
    kron_mul_elem_t!(p::Int64, q::Int64, m::Int64, a::AbstractMatrix, b::AbstractVector, c::AbstractVector)

Performs (I_p \otimes a' \otimes I_q) b, where m,n = size(a). The result is stored in c.
"""
function kron_mul_elem_t!(c::Vector, offset_c::Int64, a::AbstractMatrix, b::Vector, offset_b::Int64, p::Int64, q::Int64)
    m, n = size(a)
    length(b) >= m*p*q || throw(DimensionMismatch("The dimension of vector b, $(length(b)) doesn't correspond to order, ($p, $q)  and the dimensions of the matrix, $(size(a))"))
    length(c) >= n*p*q || throw(DimensionMismatch("The dimension of the vector c, $(length(c)) doesn't correspond to order, ($p, $q)  and the dimensions of the matrix, $(size(a))"))
    
    begin
        if p == 1 && q == 1
            # a'*b
            At_mul_B!(c, offset_c, a, 1, m, n, b, offset_b, 1)
        elseif q == 1
            #  (I_p ⊗ a')*b = vec(a'*[b_1 b_2 ... b_p])
            At_mul_B!(c, offset_c, a, 1, m, n, b, offset_b, p)
        elseif p == 1
            # (a' ⊗ I_q)*b = (b'*(a ⊗ I_q))' = vec(reshape(b,q,m)*a)
            A_mul_B!(c, offset_c, b, offset_b, q, m, a, 1, n)
        else
            # (I_p ⊗ a' ⊗ I_q)*b = vec([(a' ⊗ I_q)*b_1 (a' ⊗ I_q)*b_2 ... (a' ⊗ I_q)*b_p])
            mq = m*q
            nq = n*q
            for i=1:p
                A_mul_B!(c, offset_c, b, offset_b, q, m, a, 1, n)
                offset_b += mq
                offset_c += nq
            end
        end
    end
end

function kron_mul_elem!(c::AbstractVector, a::AbstractMatrix, b::AbstractVector, p::Int64, q::Int64)
    kron_mul_elem!(c, 1, a, b, 1, p, q)
end

#function kron_mul_elem_t!(c::Vector, a::SubArray{Float64,2,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},UnitRange{Int64}},true}, b::Vector, p::Int64, q::Int64)
#    kron_mul_elem_t!(c, , a, b, 1, p, q)
#end

function kron_mul_elem_t!(c::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true}, a::AbstractMatrix, b::SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true}, p::Int64, q::Int64)
    kron_mul_elem_t!(c.parent, c.offset1 + 1, a, b.parent, b.offset1 + 1, p, q)
end

function kron_mul_elem_t!(c::Vector, a::AbstractMatrix, b::Vector, p::Int64, q::Int64)
    kron_mul_elem_t!(c, 1, a, b, 1, p, q)
end

import Base.LinAlg: BlasInt, BlasFloat
import Base.LinAlg.BLAS: @blasfunc, libblas

function A_mul_B!(c::VecOrMat{Float64}, offset_c::Int64, a::VecOrMat{Float64},
                  offset_a::Int64, ma::Int64, na::Int64, b::VecOrMat{Float64},
                  offset_b::Int64, nb::Int64)

    gemm!('N', 'N', 1.0, Ref(a, offset_a), ma, na, Ref(b, offset_b),
          nb, 0.0, Ref(c, offset_c))
end

function A_mul_B!(c::Array{Float64,1}, offset_c::Int64, a::SubArray{Float64,2,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},UnitRange{Int64}},true}, offset_a::Int64, ma::Int64, na::Int64, b::Array{Float64,1}, offset_b::Int64, nb::Int64)
    ref_a = Ref(a, offset_a)
    ref_b = Ref(b, offset_b)
    ref_c = Ref(c, offset_c)
    lda = max(1,size(a.parent,1))
    ccall((@blasfunc(dgemm_), libblas), Void,
          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
           Ref{BlasInt}),
          'N', 'N', ma, nb,
          na, 1.0, ref_a, lda,
          ref_b, na, 0.0, ref_c,
          ma)
end

function A_mul_B!(c::Array{Float64,1}, offset_c::Int64, a::Array{Float64,1}, offset_a::Int64, ma::Int64, na::Int64, b::SubArray{Float64,2,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},UnitRange{Int64}},true}, offset_b::Int64, nb::Int64)
    ccall((@blasfunc(dgemm_), libblas), Void,
          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
           Ref{BlasInt}),
          'N', 'N', ma, nb,
          na, 1.0, Ref(a, offset_a), max(1,ma),
          Ref(b, offset_b), max(1,size(b.parent,1)), 0.0, Ref(c, offset_c),
          max(1,ma))
end

function At_mul_B!(c::VecOrMat{Float64}, offset_c::Int64, a::VecOrMat{Float64},
                  offset_a::Int64, ma::Int64, na::Int64, b::VecOrMat{Float64},
                  offset_b::Int64, nb::Int64)
    ccall((@blasfunc(dgemm_), libblas), Void,
          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
           Ref{BlasInt}),
          'T', 'N', na, nb,
          ma, 1.0, a, max(1,na),
          b, max(1,ma), 0.0, c,
          max(1,na))
end

function A_mul_Bt!(c::VecOrMat{Float64}, offset_c::Int64, a::VecOrMat{Float64},
                  offset_a::Int64, ma::Int64, na::Int64, b::VecOrMat{Float64},
                  offset_b::Int64, nb::Int64)
    gemm_t!('N', 'T', 1.0, Ref(a, offset_a), ma, na, Ref(b, offset_b),
          nb, 0.0, Ref(c, offset_c))
end

function At_mul_Bt!(c::VecOrMat{Float64}, offset_c::Int64, a::VecOrMat{Float64},
                  offset_a::Int64, ma::Int64, na::Int64, b::VecOrMat{Float64},
                  offset_b::Int64, nb::Int64)
    gemm!('T', 'T', 1.0, Ref(a, offset_a), ma, na, Ref(b, offset_b),
          nb, 0.0, Ref(c, offset_c))
end

function gemm!(ta::Char, tb::Char, alpha::Float64, a::Union{Ref{Float64},VecOrMat{Float64}},
               ma::Int64, na::Int64, b::Union{Ref{Float64},VecOrMat{Float64}}, nb::Int64,
               beta::Float64, c::Union{Ref{Float64},VecOrMat{Float64}})
    ccall((@blasfunc(dgemm_), libblas), Void,
          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
           Ref{BlasInt}),
          ta, tb, ma, nb,
          na, alpha, a, max(1,ma),
          b, max(1,na), beta, c,
          max(1,ma))
end

function gemm_t!(ta::Char, tb::Char, alpha::Float64, a::Union{Ref{Float64},VecOrMat{Float64}},
               ma::Int64, na::Int64, b::Union{Ref{Float64},VecOrMat{Float64}}, nb::Int64,
               beta::Float64, c::Union{Ref{Float64},VecOrMat{Float64}})
    ccall((@blasfunc(dgemm_), libblas), Void,
          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
           Ref{BlasInt}),
          ta, tb, ma, nb,
          na, alpha, a, max(1,ma),
          b, max(1,nb), beta, c,
          max(1,ma))
end

                      
end
