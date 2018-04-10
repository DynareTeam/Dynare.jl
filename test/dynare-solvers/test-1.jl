#=
** This script performs tests for nonlinear solvers. 
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
using TestFunctions

cd("$(rootdir)")

@testset "test unit-test-routines" begin
    @test begin
        try
            x = rosenbrock()
            f = Vector{Float64}(2)
            j = Matrix{Float64}(2,2)
            rosenbrock!(f, x)
            rosenbrock!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = powell1()
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            powell1!(f, x)
            powell1!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = powell2()
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            powell2!(f, x)
            powell2!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = wood()
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            wood!(f, x)
            wood!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = helicalvalley()
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            helicalvalley!(f, x)
            helicalvalley!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = watson(10)
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            watson!(f, x)
            watson!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = chebyquad(10)
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            chebyquad!(f, x)
            chebyquad!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = brown(10)
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            brown!(f, x)
            brown!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = discreteboundaryvalue(10)
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            discreteboundaryvalue!(f, x)
            discreteboundaryvalue!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = discreteintegralequation(10)
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            discreteintegralequation!(f, x)
            discreteintegralequation!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = trigonometric(10)
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            trigonometric!(f, x)
            trigonometric!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = variablydimensioned(10)
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            variablydimensioned!(f, x)
            variablydimensioned!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = broydentridiagonal(10)
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            broydentridiagonal!(f, x)
            broydentridiagonal!(j, x)
            true
        catch
            false
        end
    end
    @test begin
        try
            x = broydenbanded(10)
            n = length(x)
            f = Vector{Float64}(n)
            j = Matrix{Float64}(n,n)
            broydenbanded!(f, x)
            broydenbanded!(j, x)
            true
        catch
            false
        end
    end
end

cd(origdir)
