#=
** This script performs tests for the trust region nonlinear solver. 
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
using DynareSolvers
using TestFunctions

cd("$(rootdir)")

@testset "test trust-region-solver" begin
    @test begin
        x, info  = trustregion(rosenbrock!, rosenbrock!, rosenbrock(), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(powell1!, powell1!, powell1(), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(powell2!, powell2!, powell2(), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(wood!, wood!, wood(), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(helicalvalley!, helicalvalley!, helicalvalley(), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(watson!, watson!, watson(10), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(chebyquad!, chebyquad!, chebyquad(10), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(brown!, brown!, brown(10), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(discreteboundaryvalue!, discreteboundaryvalue!, discreteboundaryvalue(10), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(discreteintegralequation!, discreteintegralequation!, discreteintegralequation(10), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(trigonometric!, trigonometric!, trigonometric(10), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(variablydimensioned!, variablydimensioned!, variablydimensioned(10), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(broydentridiagonal!, broydentridiagonal!, broydentridiagonal(10), 1.0e-6, 1.0e-6, 50)
        info==1
    end
    @test begin
        x, info  = trustregion(broydenbanded!, broydenbanded!, broydenbanded(10), 1.0e-6, 1.0e-6, 50)
        info==1
    end
end

cd(origdir)
