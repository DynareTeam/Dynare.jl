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

cp("rbc1.mod", "rbc11_4.mod"; remove_destination=true)
cp("rbc1.mod", "rbc12_4.mod"; remove_destination=true)
cp("rbc1.mod", "rbc14_4.mod"; remove_destination=true)
cp("rbc2.mod", "rbc24_4.mod"; remove_destination=true)
cp("rbc1.mod", "rbc15_4.mod"; remove_destination=true)

@testset "test preprocessor-4" begin
    # Use @compile macro to compile one mod file, passing a string (name with extension).
    @testset "Test preprocessor-4-1" begin
        @test begin
            try
                @dynare "rbc11_4.mod" :savemacro :nograph
                true
            catch
                false
            end
        end
        @test isdefined(:rbc11_4)
        @test isdefined(rbc11_4, :model_) && isa(rbc11_4.model_, DynareModel.Model)
        @test isdefined(rbc11_4, :options_) && isa(rbc11_4.options_, DynareOptions.Options)
        @test isdefined(rbc11_4, :oo_) && isa(rbc11_4.oo_, DynareOutput.Output)
        @test isfile("rbc11_4.jl")
        @test isfile("rbc11_4Dynamic.jl")
        @test isfile("rbc11_4Static.jl")
        @test isfile("rbc11_4SteadyState2.jl")
        @test isfile("rbc11_4-macroexp.mod")
    end
    DynareUnitTests.clean("rbc11_4")
    rm("rbc11_4.mod")
    rm("rbc11_4-macroexp.mod")
    # Use @compile macro to compile one mod file, passing a string (name without extension).
    @testset "Test preprocessor-4-2" begin
        @test begin
            try
                @dynare "rbc12_4" :nograph :savemacro
                true
            catch
                false
            end
        end
        @test isdefined(:rbc12_4)
        @test isdefined(rbc12_4, :model_) && isa(rbc12_4.model_, DynareModel.Model)
        @test isdefined(rbc12_4, :options_) && isa(rbc12_4.options_, DynareOptions.Options)
        @test isdefined(rbc12_4, :oo_) && isa(rbc12_4.oo_, DynareOutput.Output)
        @test isfile("rbc12_4.jl")
        @test isfile("rbc12_4Dynamic.jl")
        @test isfile("rbc12_4Static.jl")
        @test isfile("rbc12_4SteadyState2.jl")
        @test isfile("rbc12_4-macroexp.mod")
    end
    DynareUnitTests.clean("rbc12_4")
    rm("rbc12_4.mod")
    rm("rbc12_4-macroexp.mod")
    # Use @compile macro to compile two mod files, passing strings (names with extensions).
    @testset "Test preprocessor-4-4" begin
        @test begin
            try
                @dynare ["rbc14_4.mod", "rbc24_4.mod"] :nograph :savemacro
                true
            catch
                false
            end
        end
        @test isdefined(:rbc14_4)
        @test isdefined(:rbc24_4)
        @test isdefined(rbc14_4, :model_) && isa(rbc14_4.model_, DynareModel.Model)
        @test isdefined(rbc14_4, :options_) && isa(rbc14_4.options_, DynareOptions.Options)
        @test isdefined(rbc14_4, :oo_) && isa(rbc14_4.oo_, DynareOutput.Output)
        @test isdefined(rbc24_4, :model_) && isa(rbc24_4.model_, DynareModel.Model)
        @test isdefined(rbc24_4, :options_) && isa(rbc24_4.options_, DynareOptions.Options)
        @test isdefined(rbc24_4, :oo_) && isa(rbc24_4.oo_, DynareOutput.Output)
        @test isfile("rbc14_4.jl")
        @test isfile("rbc24_4.jl")
        @test isfile("rbc14_4Dynamic.jl")
        @test isfile("rbc24_4Dynamic.jl")
        @test isfile("rbc14_4Static.jl")
        @test isfile("rbc24_4Static.jl")
        @test isfile("rbc14_4SteadyState2.jl")
        @test isfile("rbc14_4-macroexp.mod")
        @test isfile("rbc24_4-macroexp.mod")
    end
    DynareUnitTests.clean("rbc14_4")
    DynareUnitTests.clean("rbc24_4")
    rm("rbc14_4.mod")
    rm("rbc24_4.mod")
    rm("rbc14_4-macroexp.mod")
    rm("rbc24_4-macroexp.mod")
    # Use @compile macro to compile one mod file, passing a string (name with extension).
    @testset "Test preprocessor-4-5" begin
        @test begin
            try
                @dynare "rbc15_4.mod" :savemacro :nograph :json "compute"
                true
            catch
                false
            end
        end
        @test isdefined(:rbc15_4)
        @test isdefined(rbc15_4, :model_) && isa(rbc15_4.model_, DynareModel.Model)
        @test isdefined(rbc15_4, :options_) && isa(rbc15_4.options_, DynareOptions.Options)
        @test isdefined(rbc15_4, :oo_) && isa(rbc15_4.oo_, DynareOutput.Output)
        @test isfile("rbc15_4.jl")
        @test isfile("rbc15_4Dynamic.jl")
        @test isfile("rbc15_4Static.jl")
        @test isfile("rbc15_4SteadyState2.jl")
        @test isfile("rbc15_4-macroexp.mod")
        @test isfile("rbc15_4.json")
        @test isfile("rbc15_4_dynamic.json")
        @test isfile("rbc15_4_original.json")
        @test isfile("rbc15_4_static.json")
        @test isfile("rbc15_4_steady_state_model.json")
    end
    DynareUnitTests.clean("rbc15_4")
    rm("rbc15_4.mod")
    rm("rbc15_4-macroexp.mod")
    rm("rbc15_4.json")
    rm("rbc15_4_original.json")
    rm("rbc15_4_dynamic.json")
    rm("rbc15_4_static.json")
    rm("rbc15_4_steady_state_model.json")
end

cd(origdir)
