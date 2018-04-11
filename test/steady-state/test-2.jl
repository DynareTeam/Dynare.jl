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

cp("example2.mod", "example2_1.mod"; remove_destination=true)

@testset "test steadystate-block" begin
    @test begin
        @dynare "example2_1.mod"
        try
            steady_state!(example2_1.model_, example2_1.oo_, false)
            issteadystate(example2_1.model_, example2_1.oo_)
            true
        catch
            false
        end
    end
end

DynareUnitTests.clean("example2_1")
rm("example2_1.mod")

cd(origdir)
