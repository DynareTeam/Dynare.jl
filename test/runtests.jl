# Put Dynare's sources in the path.
if isempty(findin([abspath("../src")], LOAD_PATH))
    unshift!(LOAD_PATH, abspath("../src"))
end

# Get current directory.
rootdir0 = @__DIR__

using Base.Test

# Run tests
@testset "Dynare.jl testsuite" begin
    include("$(rootdir0)/preprocessor/test-1.jl")
    include("$(rootdir0)/preprocessor/test-2.jl")
    include("$(rootdir0)/preprocessor/test-3.jl")
    include("$(rootdir0)/preprocessor/test-4.jl")
end
