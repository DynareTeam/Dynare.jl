module KOrder

using Base.Test
#import ...taylor.FaaDiBruno: faadibruno
export make_gg!, KOrderWs

struct KOrderWs
    nstate::Integer
    nshock::Integer
    state_index::Array{Integer}
    state_range::Range
end

"""
    function make_gg!(gg,g,order,ws)

assembles the derivatives of function
gg(y,u,σ,ϵ) = [g(y,u,σ); y; u; σ; ϵ] at order 'order' 
with respect to [y, u, σ, ϵ]
"""  
function make_gg!(gg,g,order,ws)
    ngg1 = ws.nstate + 2*ws.nshock + 1
    mgg1 = ngg1 + ws.nstate
    @assert size(gg[order]) == (mgg1, ngg1^order) 
    @assert size(g[order],2) == (ws.nstate + ws.nshock + 1)^order
    @assert ws.state_range.stop <= size(g[order],1)
    if order == 1
        v2 = view(g[1],ws.state_index,:)
        copy!(gg[1],v2)
        for i = 1:(ws.nstate + 2*ws.nshock + 1)
            gg[1][ws.nstate + i,i] = 1.0
        end
    else
        n = ws.nstate + ws.nshock + 1
        i1 = CartesianIndex(1,(repmat([1], order - 1))...)
        i2 = CartesianIndex(1,(repmat([n], order - 1))...)
        i3 = ((repmat([ngg1], order))...)
        r1 = 1:n
        r2 = r1
        for i in CartesianRange(i1,i2)
            j = sub2ind(i3,(i.I)...) - 1
            v1 = view(gg[order], 1:ws.nstate, j + r1)
            v2 = view(g[order], ws.state_index, r2)
            copy!(v1,v2)
            r2 += n
        end
    end
end

end
