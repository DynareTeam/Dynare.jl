#=
** This script performs tests for the steady state numerical solver. 
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
using SteadyState
using DynareUnitTests

cd("$(rootdir)")

cp("example1.mod", "example1_1.mod"; remove_destination=true)

steadystate = [1.080683, 0.803593, 11.083609, 0.000000, 0.291756, 0.000000]

@testset "test steadystate-solver" begin
    @test begin
        try
            @dynare "example1_1.mod"
            true
        catch
            false
        end
    end
    model = example1_1.model_
    oo = example1_1.oo_
    @test begin
        try
            yinit = [1.0, .8, 10.0, 0.0, .3, 0.1]
            steady_state(model, oo, yinit, false)
            true
        catch
            false
        end
    end
    @test begin
        try
            yinit = [2.0, .8, 6.0, 0.0, .3, 0.1]
            steady!(model, oo, yinit)
            all((oo.steady_state-steadystate).<1.0e-6)
        catch
            false
        end
    end
    @test begin
        try
            yinit = [.5, .3, 12.0, 0.0, .2, -0.1]
            steady!(model, oo, yinit)
            all((oo.steady_state-steadystate).<1.0e-6)
        catch
            false
        end
    end
    @test begin
        try
            yinit = [.2, .3, 4.0, 0.0, .2, -0.1]
            steady!(model, oo, yinit)
            all((oo.steady_state-steadystate).<1.0e-6)
        catch
            false
        end
    end
    @test begin
        try
            yinit = [.1, .3, 4.0, 0.0, .3, 0.1]
            steady!(model, oo, yinit)
            all((oo.steady_state-steadystate).<1.0e-6)
        catch
            false
        end
    end
    @test begin
        try
            yinit = [20.0, .1, 3.0, 0.0, .3, 0.1]
            steady!(model, oo, yinit)
            all((oo.steady_state-steadystate).<1.0e-6)
        catch
            false
        end
    end
end

DynareUnitTests.clean("example1_1")
rm("example1_1.mod")

cd(origdir)
