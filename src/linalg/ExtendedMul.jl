module ExtendedMul

import Base.LinAlg: A_mul_B!, At_mul_B!, A_mul_Bt!, At_mul_Bt
import Base.BLAS: gemm!

import Base.LinAlg: BlasInt, BlasFloat
import Base.LinAlg.BLAS: @blasfunc, libblas

export A_mul_B!, At_mul_B!, A_mul_Bt!, At_mul_B!

function A_mul_B!(c::VecOrMat{Float64}, offset_c::Int64, a::VecOrMat{Float64},
                  offset_a::Int64, ma::Int64, na::Int64, b::VecOrMat{Float64},
                  offset_b::Int64, nb::Int64)

    gemm!('N', 'N', 1.0, Ref(a, offset_a), ma, na, Ref(b, offset_b),
          nb, 0.0, Ref(c, offset_c))
end


function A_mul_B!(c::Array{Float64,1}, offset_c::Int64, a::SubArray{Float64,2,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},UnitRange{Int64}},true}, offset_a::Int64, ma::Int64, na::Int64, b::Array{Float64,1}, offset_b::Int64, nb::Int64)
    if offset_a != 1
        throw(DimensionMismatch("offset_a must be 1"))
    end
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
    if offset_b != 1
        throw(DimensionMismatch("offset_a must be 1"))
    end
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
          'T', 'N', ma, nb,
          na, 1.0, Ref(a, offset_a), max(1,na),
          Ref(b, offset_b), max(1,na), 0.0, Ref(c, offset_c),
          max(1,ma))
end

function A_mul_Bt!(c::VecOrMat{Float64}, offset_c::Int64, a::VecOrMat{Float64},
                  offset_a::Int64, ma::Int64, na::Int64, b::VecOrMat{Float64},
                  offset_b::Int64, nb::Int64)
    ccall((@blasfunc(dgemm_), libblas), Void,
          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
           Ref{BlasInt}),
          'N', 'T', ma, nb,
          na, 1.0, Ref(a, offset_a), max(1,ma),
          Ref(b, offset_b), max(1,nb), 0.0, Ref(c, offset_c),
          max(1,ma))
end

function At_mul_Bt!(c::VecOrMat{Float64}, offset_c::Int64, a::VecOrMat{Float64},
                  offset_a::Int64, ma::Int64, na::Int64, b::VecOrMat{Float64},
                  offset_b::Int64, nb::Int64)

    ccall((@blasfunc(dgemm_), libblas), Void,
          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
           Ref{BlasInt}),
          'T', 'T', ma, nb,
          na, 1.0, Ref(a, offset_a), max(1,na),
          Ref(b, offset_b), max(1,nb), 0.0, Ref(c, offset_c),
          max(1,ma))
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
end

                      
end
