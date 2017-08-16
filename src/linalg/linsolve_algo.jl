module linsolve_algo

import Base.LinAlg.BlasInt
import Base.LinAlg.BLAS.@blasfunc
import Base.LinAlg.BLAS.libblas
import Base.LinAlg.LAPACK: liblapack, chklapackerror

export LinSolveWS, linsolve_core!

struct LinSolveWS
    ipiv::Vector{BlasInt}

    function LinSolveWS(n)
        ipiv = Vector{BlasInt}(n)
        new(ipiv)
    end
end

function linsolve_core!(ws::LinSolveWS,trans::Ref{UInt8},a::StridedMatrix{Float64},b::StridedVecOrMat{Float64})
    mm,nn = size(a)
    m = Ref{BlasInt}(mm)
    n = Ref{BlasInt}(nn)
    nhrs = Ref{BlasInt}(size(b,2))
    lda = Ref{BlasInt}(max(1,stride(a,2)))
    ldb = Ref{BlasInt}(max(1,stride(b,2)))
    info = Ref{BlasInt}(0)

    lu!(a,ws.ipiv)
    ccall((@blasfunc(dgetrs_), liblapack), Void,
          (Ref{UInt8},Ref{BlasInt},Ref{BlasInt},Ptr{Float64},Ref{BlasInt},
           Ptr{BlasInt},Ptr{Float64},Ref{BlasInt},Ref{BlasInt}),
          trans,n,nhrs,a,lda,ws.ipiv,b,ldb,info)
    if info[] != 0
        println("dgetrs ",info[])
        chklapackerror(info[])
    end
end

function linsolve_core!(ws::LinSolveWS,trans::Ref{UInt8},a::StridedMatrix{Float64},b::StridedVecOrMat{Float64},c::StridedVecOrMat{Float64})
    mm,nn = size(a)
    m = Ref{BlasInt}(mm)
    n = Ref{BlasInt}(nn)
    lda = Ref{BlasInt}(max(1,stride(a,2)))
    ldb = Ref{BlasInt}(max(1,stride(b,2)))
    ldc = Ref{BlasInt}(max(1,stride(c,2)))
    info = Ref{BlasInt}(0)

    lu!(a,ws.ipiv)
    nhrs = Ref{BlasInt}(size(b,2))
    ccall((@blasfunc(dgetrs_), liblapack), Void,
          (Ref{UInt8},Ref{BlasInt},Ref{BlasInt},Ptr{Float64},Ref{BlasInt},
           Ptr{BlasInt},Ptr{Float64},Ref{BlasInt},Ref{BlasInt}),
          trans,n,nhrs,a,lda,ws.ipiv,b,ldb,info)
    if info[] != 0
        println("dgetrs ",info[])
        chklapackerror(info[])
    end
    nhrs = Ref{BlasInt}(size(c,2))
    ccall((@blasfunc(dgetrs_), liblapack), Void,
          (Ref{UInt8},Ref{BlasInt},Ref{BlasInt},Ptr{Float64},Ref{BlasInt},
           Ptr{BlasInt},Ptr{Float64},Ref{BlasInt},Ref{BlasInt}),
          trans,n,nhrs,a,lda,ws.ipiv,c,ldc,info)
    if info[] != 0
        println("dgetrs ",info[])
        chklapackerror(info[])
    end
end

function linsolve_core!(ws::LinSolveWS,trans::Ref{UInt8},a::StridedMatrix{Float64},b::StridedVecOrMat{Float64},c::StridedVecOrMat{Float64},d::StridedVecOrMat{Float64})
    mm,nn = size(a)
    m = Ref{BlasInt}(mm)
    n = Ref{BlasInt}(nn)
    nhrs = Ref{BlasInt}(size(b,2))
    lda = Ref{BlasInt}(max(1,stride(a,2)))
    ldb = Ref{BlasInt}(max(1,stride(b,2)))
    ldc = Ref{BlasInt}(max(1,stride(c,2)))
    ldd = Ref{BlasInt}(max(1,stride(d,2)))
    info = Ref{BlasInt}(0)

    lu!(a,ws.ipiv)
    ccall((@blasfunc(dgetrs_), liblapack), Void,
          (Ref{UInt8},Ref{BlasInt},Ref{BlasInt},Ptr{Float64},Ref{BlasInt},
           Ptr{BlasInt},Ptr{Float64},Ref{BlasInt},Ref{BlasInt}),
          trans,n,nhrs,a,lda,ws.ipiv,b,ldb,info)
    if info[] != 0
        println("dgetrs ",info[])
        chklapackerror(info[])
    end
    ccall((@blasfunc(dgetrs_), liblapack), Void,
          (Ref{UInt8},Ref{BlasInt},Ref{BlasInt},Ptr{Float64},Ref{BlasInt},
           Ptr{BlasInt},Ptr{Float64},Ref{BlasInt},Ref{BlasInt}),
          trans,n,nhrs,a,lda,ws.ipiv,c,ldc,info)
    if info[] != 0
        println("dgetrs ",info[])
        chklapackerror(info[])
    end
    ccall((@blasfunc(dgetrs_), liblapack), Void,
          (Ref{UInt8},Ref{BlasInt},Ref{BlasInt},Ptr{Float64},Ref{BlasInt},
           Ptr{BlasInt},Ptr{Float64},Ref{BlasInt},Ref{BlasInt}),
          trans,n,nhrs,a,lda,ws.ipiv,d,ldd,info)
    if info[] != 0
        println("dgetrs ",info[])
        chklapackerror(info[])
    end
end

function lu!(a,ipiv)
    mm,nn = size(a)
    m = Ref{BlasInt}(mm)
    n = Ref{BlasInt}(nn)
    lda = Ref{BlasInt}(max(1,stride(a,2)))
    info = Ref{BlasInt}(0)
    ccall((@blasfunc(dgetrf_), liblapack), Void,
          (Ref{BlasInt},Ref{BlasInt},Ptr{Float64},Ref{BlasInt},
           Ptr{BlasInt},Ref{BlasInt}),
          m,n,a,lda,ipiv,info)
    if info[] != 0
        println("dgetrf ",info[])
        chklapackerror(info[])
    end
end

end
