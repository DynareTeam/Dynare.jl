if isempty(findin([abspath("../src")], LOAD_PATH))
    unshift!(LOAD_PATH, abspath("../src"))
end

using Dynare

include("ramst.jl")
include("example1.jl")
