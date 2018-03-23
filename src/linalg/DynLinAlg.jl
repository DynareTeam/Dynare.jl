module DynLinAlg

include("KroneckerUtils.jl")
using .KroneckerUtils
include("LinSolveAlgo.jl")
using .LinSolveAlgo
include("QrAlgo.jl")
using .QrAlgo
# QuasiUpperTriangular Matrix of real numbers
include("QUT.jl")
using .QUT
include("SchurAlgo.jl")
using .SchurAlgo 

end

