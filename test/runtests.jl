# Put Dynare's sources in the path.
if isempty(findin([abspath("../src")], LOAD_PATH))
    unshift!(LOAD_PATH, abspath("../src"))
end

# Get current directory.
rootdir = @__DIR__

using Base.Test

# Run tests
@testset "Dynare.jl testsuite" begin
    include("$(rootdir)/preprocessor/test-1.jl")
end
