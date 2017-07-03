module LinAlg

# Kronecker operations
include("kronecker_utils.jl")

# linear solver
include("linsolve_algo.jl")

# QR 
include("qr_algo.jl")

# QuasiUpperTriangular type
include("quasi_upper_triangular.jl")

# Schur decompositions
include("schur.jl")
export DgeesWS, dgees!, DggesWS, dgges!

end

