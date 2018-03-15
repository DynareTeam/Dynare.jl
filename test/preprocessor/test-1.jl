#=
** This script performs tests for the @compile macro. 
**
** Copyright (C) 2018 Dynare Team
**
** This file is part of Dynare.
**
** Dynare is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** Dynare is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
=#

    
# Put Dynare's sources in the path
if isempty(findin([abspath("../../src")], LOAD_PATH))
    unshift!(LOAD_PATH, abspath("../../src"))
end

using Base.Test
using Dynare
using DynareUnitTests

rootdir = @__DIR__
origdir = pwd()

cd("$(rootdir)")

cp("rbc1.mod", "rbc11.mod"; remove_destination=true)
cp("rbc1.mod", "rbc12.mod"; remove_destination=true)
cp("rbc1.mod", "rbc13.mod"; remove_destination=true)
cp("rbc1.mod", "rbc15.mod"; remove_destination=true)
cp("rbc2.mod", "rbc25.mod"; remove_destination=true)
cp("rbc1.mod", "rbc16.mod"; remove_destination=true)
cp("rbc2.mod", "rbc26.mod"; remove_destination=true)

@testset "test preprocessor-1" begin
    # Use @compile macro to compile one mod file, passing a string (name with extension).
    @testset "Test preprocessor-1-1" begin
        @test @compile "rbc11.mod"
        @test Base.isfile("rbc11.jl")
        @test isfile("rbc11Dynamic.jl")
        @test isfile("rbc11Static.jl")
        @test isfile("rbc11SteadyState2.jl")
    end
    DynareUnitTests.clean("rbc11")
    rm("rbc11.mod")
    # Use @compile macro to compile one mod file, passing a string (name without extension).
    @testset "Test preprocessor-1-2" begin
        @test @compile "rbc12"
        @test isfile("rbc12.jl")
        @test isfile("rbc12Dynamic.jl")
        @test isfile("rbc12Static.jl")
        @test isfile("rbc12SteadyState2.jl")
    end
    DynareUnitTests.clean("rbc12")
    rm("rbc12.mod")
    # Use @compile macro to compile one mod file, passing a symbol.
    @testset "Test preprocessor-1-3" begin
        @test @compile :rbc13
        @test isfile("rbc13.jl")
        @test isfile("rbc13Dynamic.jl")
        @test isfile("rbc13Static.jl")
        @test isfile("rbc13SteadyState2.jl")
    end
    DynareUnitTests.clean("rbc13")
    rm("rbc13.mod")
    # Use compile function to compile one mod file, passing a string (name with a wrong extension, this should fail).
    @testset "Test preprocessor-1-4" begin
        @test_throws Core.ErrorException Dynare.compile("rbc3.zut")
        @test ~isfile("rbc3.jl")
        @test ~isfile("rbc3Dynamic.jl")
        @test ~isfile("rbc3Static.jl")
        @test ~isfile("rbc3SteadyState2.jl")
    end
    # Use @compile macro to compile two mod files, passing strings (names with extensions).
    @testset "Test preprocessor-1-5" begin
        @test @compile "rbc15.mod" "rbc25.mod"
        @test isfile("rbc15.jl")
        @test isfile("rbc25.jl")
        @test isfile("rbc15Dynamic.jl")
        @test isfile("rbc25Dynamic.jl")
        @test isfile("rbc15Static.jl")
        @test isfile("rbc25Static.jl")
        @test isfile("rbc15SteadyState2.jl")
    end
    DynareUnitTests.clean("rbc15")
    DynareUnitTests.clean("rbc25")
    rm("rbc15.mod")
    rm("rbc25.mod")
    # Use @compile macro to compile two mod files, passing symbols.
    @testset "Test preprocessor-1-6" begin
        @test @compile :rbc16 :rbc26
        @test isfile("rbc16.jl")
        @test isfile("rbc26.jl")
        @test isfile("rbc16Dynamic.jl")
        @test isfile("rbc26Dynamic.jl")
        @test isfile("rbc16Static.jl")
        @test isfile("rbc26Static.jl")
        @test isfile("rbc16SteadyState2.jl")
    end
    DynareUnitTests.clean("rbc16")
    DynareUnitTests.clean("rbc26")
    rm("rbc16.mod")
    rm("rbc26.mod")
end

cd(origdir)
