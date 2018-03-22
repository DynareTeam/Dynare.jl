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

rootdir = @__DIR__
origdir = pwd()

if isempty(findin([abspath("$(rootdir)/../../src")], LOAD_PATH))
    unshift!(LOAD_PATH, abspath("$(rootdir)/../../src"))
end

using Base.Test
using Dynare
using DynareUnitTests

cd("$(rootdir)")

cp("rbc1.mod", "rbc11.mod"; remove_destination=true)
cp("rbc1.mod", "rbc12.mod"; remove_destination=true)
cp("rbc1.mod", "rbc14.mod"; remove_destination=true)
cp("rbc2.mod", "rbc24.mod"; remove_destination=true)

@testset "test preprocessor-1" begin
    # Use @compile macro to compile one mod file, passing a string (name with extension).
    @testset "Test preprocessor-1-1" begin
        @test begin
            try
                @compile "rbc11.mod"
                true
            catch
                false
            end
        end
        @test Base.isfile("rbc11.jl")
        @test isfile("rbc11Dynamic.jl")
        @test isfile("rbc11Static.jl")
        @test isfile("rbc11SteadyState2.jl")
    end
    DynareUnitTests.clean("rbc11")
    rm("rbc11.mod")
    # Use @compile macro to compile one mod file, passing a string (name without extension).
    @testset "Test preprocessor-1-2" begin
        @test begin
            try
                @compile "rbc12"
                true
            catch
                false
            end
        end
        @test isfile("rbc12.jl")
        @test isfile("rbc12Dynamic.jl")
        @test isfile("rbc12Static.jl")
        @test isfile("rbc12SteadyState2.jl")
    end
    DynareUnitTests.clean("rbc12")
    rm("rbc12.mod")
    # Use compile function to compile one mod file, passing a string (name with a wrong extension, this should fail).
    @testset "Test preprocessor-1-3" begin
        @test begin
            try
                Dynare.compile("rbc3.zut")
                false
            catch
                true
            end
        end
        @test ~isfile("rbc3.jl")
        @test ~isfile("rbc3Dynamic.jl")
        @test ~isfile("rbc3Static.jl")
        @test ~isfile("rbc3SteadyState2.jl")
    end
    # Use @compile macro to compile two mod files, passing strings (names with extensions).
    @testset "Test preprocessor-1-4" begin
        @test begin
            try
                @compile ["rbc14.mod", "rbc24.mod"]
                true
            catch
                false
            end
        end
        @test isfile("rbc14.jl")
        @test isfile("rbc24.jl")
        @test isfile("rbc14Dynamic.jl")
        @test isfile("rbc24Dynamic.jl")
        @test isfile("rbc14Static.jl")
        @test isfile("rbc24Static.jl")
        @test isfile("rbc14SteadyState2.jl")
    end
    DynareUnitTests.clean("rbc14")
    DynareUnitTests.clean("rbc24")
    rm("rbc14.mod")
    rm("rbc24.mod")
end

cd(origdir)
