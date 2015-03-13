import Base.LinAlg.BlasInt
import Base.blasfunc

# Equivalent of mjdgges
# S and T are modified in place
@eval begin
    function qz!(S::Matrix{Float64}, T::Matrix{Float64})
        function selctg(alphar::Ptr{Float64}, alphai::Ptr{Float64}, beta::Ptr{Float64})
            const criterium::Float64 = 1.0 + 1e-6
            ai = unsafe_load(alphai, 1)
            ar = unsafe_load(alphar, 1)
            b = unsafe_load(beta, 1)
            
            convert(BlasInt, ar^2 + ai^2 < criterium * b^2)
        end
        
        selctg2 = cfunction(selctg, BlasInt, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}))
        (n,) = size(S)
        @assert size(S) == (n,n)
        @assert size(T) == (n,n)
        @assert stride(S, 1) == 1
        @assert stride(T, 1) == 1
        
        alphar = Array(Float64, n)
        alphai = Array(Float64, n)
        beta = Array(Float64, n)
        lwork = 16*n+16 # Same heuristic choice than mjdgges
        work = Array(Float64, lwork)
        vsl = Array(Float64, n, n)
        vsr = Array(Float64, n, n)
        bwork = Array(BlasInt, n)
        sdim = Array(BlasInt, 1)
        info = Array(BlasInt, 1)

        ccall( ( $(blasfunc(:dgges_)), Base.liblapack_name), Void, (Ptr{Uint8}, Ptr{Uint8}, Ptr{Uint8}, Ptr{Void},
                                                                    Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                                                                    Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
                                                                    Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                                                                    Ptr{Float64}, Ptr{BlasInt},
                                                                    Ptr{Float64}, Ptr{BlasInt},
                                                                    Ptr{Float64}, Ptr{BlasInt},
                                                                    Ptr{BlasInt}, Ptr{BlasInt}),
              &'N', &'V', &'S', selctg2, &n, S, &stride(S, 2), T, &stride(T, 2),
              sdim, alphar, alphai, beta, vsl, &stride(vsl, 2),
              vsr, &stride(vsr, 2), work, &lwork, bwork, info)
        
        @assert info[1] == 0
        
        eigs = Array(Complex128, n)
        for i = 1:n
            eigs[i] = alphar[i]/beta[i] + (alphai[i] == 0 && beta[i] == 0 ? 0.0 : alphai[i]/beta[i])*im
        end
        
        (vsr, sdim[1], eigs)
    end
end
