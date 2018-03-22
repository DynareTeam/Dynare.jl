#=
** This script performs tests for the @compile macro (with options). 
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

rootdir = @__DIR__
origdir = pwd()

if isempty(findin([abspath("$(rootdir)/../../src")], LOAD_PATH))
    unshift!(LOAD_PATH, abspath("$(rootdir)/../../src"))
end

using Base.Test
using Dynare
using DynareUnitTests

cd("$(rootdir)")

cp("rbc1.mod", "rbc11_3.mod"; remove_destination=true)
cp("rbc1.mod", "rbc12_3.mod"; remove_destination=true)
cp("rbc1.mod", "rbc14_3.mod"; remove_destination=true)
cp("rbc2.mod", "rbc24_3.mod"; remove_destination=true)

@testset "test preprocessor-3" begin
    # Use @compile macro to compile one mod file, passing a string (name with extension).
    @testset "Test preprocessor-3-1" begin
        @test begin
            try
                @compile "rbc11_3.mod" @PreprocessorOptions :savemacro :nograph
                true
            catch
                false
            end
        end
        @test isfile("rbc11_3.jl")
        @test isfile("rbc11_3Dynamic.jl")
        @test isfile("rbc11_3Static.jl")
        @test isfile("rbc11_3SteadyState2.jl")
        @test isfile("rbc11_3-macroexp.mod")
    end
    DynareUnitTests.clean("rbc11_3")
    rm("rbc11_3.mod")
    rm("rbc11_3-macroexp.mod")
    # Use @compile macro to compile one mod file, passing a string (name without extension).
    @testset "Test preprocessor-3-2" begin
        @test begin
            try
                @compile "rbc12_3" @PreprocessorOptions :nograph :savemacro
                true
            catch
                false
            end
        end
        @test isfile("rbc12_3.jl")
        @test isfile("rbc12_3Dynamic.jl")
        @test isfile("rbc12_3Static.jl")
        @test isfile("rbc12_3SteadyState2.jl")
        @test isfile("rbc12_3-macroexp.mod")
    end
    DynareUnitTests.clean("rbc12_3")
    rm("rbc12_3.mod")
    rm("rbc12_3-macroexp.mod")
    # Use @compile macro to compile two mod files, passing strings (names with extensions).
    @testset "Test preprocessor-1-4" begin
        @test begin
            try
                @compile ["rbc14_3.mod", "rbc24_3.mod"] @PreprocessorOptions :nograph :savemacro
                true
            catch
                false
            end
        end
        @test isfile("rbc14_3.jl")
        @test isfile("rbc24_3.jl")
        @test isfile("rbc14_3Dynamic.jl")
        @test isfile("rbc24_3Dynamic.jl")
        @test isfile("rbc14_3Static.jl")
        @test isfile("rbc24_3Static.jl")
        @test isfile("rbc14_3SteadyState2.jl")
        @test isfile("rbc14_3-macroexp.mod")
        @test isfile("rbc24_3-macroexp.mod")
    end
    DynareUnitTests.clean("rbc14_3")
    DynareUnitTests.clean("rbc24_3")
    rm("rbc14_3.mod")
    rm("rbc24_3.mod")
    rm("rbc14_3-macroexp.mod")
    rm("rbc24_3-macroexp.mod")
end

cd(origdir)
