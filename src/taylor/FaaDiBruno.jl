module FaaDiBruno
using Combinatorics
import Base: start, next, done
import KroneckerUtils: a_mul_kron_b!

export FaaDiBrunoWs, faa_di_bruno!, partial_faa_di_bruno!

ttuple = Tuple{ Array{Int,1}, Array{Array{Int,1},1} }
tatuple = Array{ttuple,1}
trecipees = Array{Array{tatuple,1},1}

struct FaaDiBrunoWs
    recipees::Array{Array{tatuple,1},1}
    work1::Vector{Float64}
    work2::Vector{Float64}
    matrixvector::Array{AbstractArray,1}
    function FaaDiBrunoWs(neq,nvar,recipees,order)
        setup_recipees!(recipees,order)
        work1 = Vector{Float64}(neq*nvar^order)
        work2 = Vector{Float64}(neq*nvar^order)
        new(recipees,work1,work2)
    end
end

function FaaDiBrunoWs(neq,nvar,order)
    recipees = Array{Array{tatuple,1},1}(order)
    FaaDiBrunoWs(neq,nvar,recipees,order)
end

multicat(a::Array{Array{Int64,1},1}) =  [a[i][j] for i in eachindex(a) for j in eachindex(a[i])]
                  
function setup_recipees!(recipees::trecipees,order::Int)
    for i in 1:order
        r1 = Array{tatuple,1}(0)
        p1 = integer_partitions(i)
        for j in 1:i
            r2 = Array{tatuple,1}(0)
            c = tatuple(0)
            p2 = collect(partitions(collect(1:i),j))
            map!(x -> sort(x,by=length),p2,p2)
            patterns = map(x->map(y->length(y),x),p2)
            p1j = filter(x->length(x) == j,p1)
            for k = 1:length(p1j)
                p1jk = p1j[k]
                r2 = Array{Array{Int,1},1}(0)
                for m = 1:length(p2)
                    if patterns[m] == sort(p1jk)
                        v = multicat(p2[m])
                        push!(r2,ipermute!(collect(1:length(v)),v))
                    end
                end
                push!(c,(p1jk,r2))
            end
            push!(r1,c)
        end
        recipees[i] = r1
    end
end

function faa_di_bruno!(dfg::AbstractMatrix{Float64},f::Array{Matrix{Float64},1},g::Array{Matrix{Float64}},order::Int,ws::FaaDiBrunoWs)
    m = size(f[1],1)
    n = size(g[1],2)
    work1 = reshape(view(ws.work1,1:m*n^order),m,n^order)
    if order == 1
        A_mul_B!(dfg,f[1],g[1])
    elseif order == 2
        A_mul_B!(dfg,f[1],g[2])
        a_mul_kron_b!(work1,f[2],g[1],2)
        dfg .+= work1
    else
        for i = 1:order
            if i == 1 && g[order] != [0]
                A_mul_B!(dfg,f[1],g[order])
#            elseif i == order
#                a_mul_kron_b!(work1,f[i],g[1],i)
#                dfg .+= work1
#                println(work1[1])
#                println(dfg[1])
            else
                apply_recipees!(dfg,ws.recipees[order][i],f[i],g,order,ws)
            end
        end
    end
end

"""
    function partial_faa_di_bruno!(dfg,f,g,order,ws)
computes the derivatives of f(g()) at order "order"
but without the term involving order'th derivative of g
"""
function partial_faa_di_bruno!(dfg::AbstractMatrix{Float64},f::Array{Matrix{Float64},1},g::Array{Matrix{Float64}},order::Int,ws::FaaDiBrunoWs)
    m = size(f[1],1)
    n = size(g[1],2)
    if order == 1
        throw(ArgumentError("Can't run partial_faa_di_bruno() for order == 1"))
    elseif order == 2
        a_mul_kron_b!(dfg,f[2],g[1],2,ws.work1,ws.work2)
    else
        for i = 1:order
            apply_recipees!(dfg,ws.recipees[order][i],f[i],g,order,ws)
        end
    end
end

function apply_recipees!(dfg::Matrix{Float64},recipees::tatuple,f::AbstractArray,g::Array{Array{Float64,2},1},order::Int64,ws::FaaDiBrunoWs)
    m = size(f,1)
    n = size(g[1],2)
    dims = (m, repmat([n],order)...)
    work1 = reshape(view(ws.work1,1:m*n^order),m,n^order)
    for i = 1:length(recipees)
        recipees1 = recipees[i][1]
        recipees2 = recipees[i][2]
        fill!(work1,0.0)
        fill!(ws.work2,0.0)
        a_mul_kron_b!(work1,f,g[recipees1],ws.work2)
        a = reshape(work1, dims)
        d = reshape(dfg, dims)
        p = length(recipees2)
        matrixvector = Array{AbstractArray}(p)
        dims0 = (1,recipees2[1]+1...)
        for j = 1:length(recipees2)
            dims1 = (1,recipees2[j]+1...)
            matrixvector[j] = PermutedDimsArray(a,dims1)
        end
        add!(d,matrixvector)
    end
end

"""
    function add!(b::Matrix{Float64},a::Array{AbstractArray,1})
adds to matrix b the sum of the (permuted) matrices in vector of matrices a

The inside loop is unfolded up to 15 matrices
"""
function add!(b::AbstractArray{Float64},a::Array{AbstractArray,1})
    n = length(a)
    @inbounds b .+= a[1]
    if n == 1
        return
    elseif n == 2
        @inbounds b .+= a[2]
    elseif n == 3
        a2 = a[2]
        a3 = a[3]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i]
        end
    elseif n == 4
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i]
        end
    elseif n == 5
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i]
        end
    elseif n == 6
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        a6 = a[6]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i] + a6[i]
        end
    elseif n == 7
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        a6 = a[6]
        a7 = a[7]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i]
        end
    elseif n == 8
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        a6 = a[6]
        a7 = a[7]
        a8 = a[8]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] + a8[i]
        end
    elseif n == 9
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        a6 = a[6]
        a7 = a[7]
        a8 = a[8]
        a9 = a[9]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] + a8[i] + a9[i]
        end
    elseif n == 10
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        a6 = a[6]
        a7 = a[7]
        a8 = a[8]
        a9 = a[9]
        a10 = a[10]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] + a8[i] + a9[i] + a10[i]
        end
    elseif n == 11
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        a6 = a[6]
        a7 = a[7]
        a8 = a[8]
        a9 = a[9]
        a10 = a[10]
        a11 = a[11]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] + a8[i] + a9[i] + a10[i] + a11[i]
        end
    elseif n == 12
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        a6 = a[6]
        a7 = a[7]
        a8 = a[8]
        a9 = a[9]
        a10 = a[10]
        a11 = a[11]
        a12 = a[12]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] + a8[i] + a9[i] + a10[i] + a11[i] + a12[i]
        end
    elseif n == 13
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        a6 = a[6]
        a7 = a[7]
        a8 = a[8]
        a9 = a[9]
        a10 = a[10]
        a11 = a[11]
        a12 = a[12]
        a13 = a[13]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] + a8[i] + a9[i] + a10[i] + a11[i] + a12[i] + a13[i]
        end
    elseif n == 14
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        a6 = a[6]
        a7 = a[7]
        a8 = a[8]
        a9 = a[9]
        a10 = a[10]
        a11 = a[11]
        a12 = a[12]
        a13 = a[13]
        a14 = a[14]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] + a8[i] + a9[i] + a10[i] + a11[i] + a12[i] + a13[i] + a14[i]
        end
    elseif n == 15
        a2 = a[2]
        a3 = a[3]
        a4 = a[4]
        a5 = a[5]
        a6 = a[6]
        a7 = a[7]
        a8 = a[8]
        a9 = a[9]
        a10 = a[10]
        a11 = a[11]
        a12 = a[12]
        a13 = a[13]
        a14 = a[14]
        a15 = a[15]
        @simd for i in CartesianRange(size(b))
            @inbounds b[i] += a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] + a8[i] + a9[i] + a10[i] + a11[i] + a12[i] + a13[i] + a14[i] + a15[i]
        end
    else
        na = length(a)
        @inbounds for i in CartesianRange(size(b))
            @simd for j = 2:na
                b[i] += a[j][i]
            end
        end
    end
end
              
        
end

