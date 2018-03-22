#=
** This script performs tests for the PreprocessorOptions type. 
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

@testset "test preprocessor-2" begin
    @test begin
        # Test that we can instantiate an object.
        try
            opt  = Dynare.PreprocessorOptions()
            true
        catch
            false
        end
    end
    @test begin
        # Test that we cannot access a field that does not exist.
        try
            opt  = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :noddy, "toyland")
            false
        catch
            true
        end
    end
    @test begin
        # Test that we cannot assign a non expected type to a member (field).
        try
            opt  = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :debug, "true")
            false
        catch
            true
        end
    end
    @test begin
        # Test that we can change the value of an option.
        try
            opt  = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :debug, true)
            opt.debug==true
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (default option values).
        try
            opt = Dynare.PreprocessorOptions()
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed debug option).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :debug, true)
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "nopreprocessoroutput") && contains(str, "debug")
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed debug and noclearall options).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :debug, true)
            Dynare.PreprocessorOptions!(opt, :noclearall, true)
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "nopreprocessoroutput") && contains(str, "debug") && contains(str, "noclearall")
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed option savemacro).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :savemacro, "toto.mod")
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "nopreprocessoroutput") && contains(str, "savemacro=toto.mod")
        catch
            false
        end
    end
    for option in [:savemacro :onlymacro :nolinemacro :noemptylinemacro :notmpterms :nolog :warn_uninit :nograph :nointeractive :parallel :parallel_slave_open_mode :parallel_test :nostrict :stochastic :fast :minimal_workspace :compute_xrefs]
        @test begin
            # Test that we can print the list of options in a string.
            try
                opt = Dynare.PreprocessorOptions()
                Dynare.PreprocessorOptions!(opt, option, true)
                str = Dynare.print(opt)
                contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "nopreprocessoroutput") && contains(str, "$(string(option))")
            catch
                false
            end
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed option parallel).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :parallel, "mycluster")
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "nopreprocessoroutput") && contains(str, "parallel=mycluster")
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed option conffile).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :conffile, "/home/toto/.dynare.m")
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "nopreprocessoroutput") && contains(str, "conffile=/home/toto/.dynare.m")
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed option conffile).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :conffile, "/home/toto/.dynare.m")
            str = Dynare.print(opt)
            str==" conffile=/home/toto/.dynare.m output=dynamic language=julia nopreprocessoroutput" 
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed option output).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :output, "second")
            str = Dynare.print(opt)
            contains(str, "output=second") && contains(str, "language=julia") && contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed option language).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :language, "c++")
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && contains(str, "language=c++") && contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed option language).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :language, "matlab")
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && !contains(str, "language") && contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed option param_derivs_order).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :params_derivs_order, 1)
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "params_derivs_order=1") && contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
    @test begin
        # Test that we can print the list of options in a string (provide wrong type to option param_derivs_order).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :params_derivs_order, "one")
            false
        catch
            true
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed option json).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :json, "transform")
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "json=transform") && contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
    for option in [:jsonstdout :onlyjson :jsonderivsimple :nopathchange]
        @test begin
            # Test that we can print the list of options in a string.
            try
                opt = Dynare.PreprocessorOptions()
                Dynare.PreprocessorOptions!(opt, option, true)
                str = Dynare.print(opt)
                contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "$(string(option))") && contains(str, "nopreprocessoroutput")
            catch
                false
            end
        end
    end
    @test begin
        # Test that we can print the list of options in a string (changed nopreprocessoroutput option).
        try
            opt = Dynare.PreprocessorOptions()
            Dynare.PreprocessorOptions!(opt, :nopreprocessoroutput, false)
            str = Dynare.print(opt)
            contains(str, "output=dynamic") && contains(str, "language=julia") && !contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
    @test begin
        # Test @PreprocessorOptions (added one option)
        try
            options = Dynare.@PreprocessorOptions :debug
            str = Dynare.print(options)
            contains(str, "debug") && contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
    @test begin
        # Test @PreprocessorOptions (added two options)
        try
            options = Dynare.@PreprocessorOptions :debug :nograph
            str = Dynare.print(options)
            contains(str, "debug") && contains(str, "nograph") && contains(str, "output=dynamic") && contains(str, "language=julia") && contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
    @test begin
        # Test @PreprocessorOptions (added two options and changed lanquage value)
        try
            options = Dynare.@PreprocessorOptions :debug :nograph :language "c++"
            str = Dynare.print(options)
            contains(str, "debug") && contains(str, "nograph") && contains(str, "output=dynamic") && contains(str, "language=c++") && contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
    @test begin
        # Test @PreprocessorOptions (added two options, changed lanquage value and removed nopreprocessoroutput)
        try
            options = Dynare.@PreprocessorOptions :debug :nograph :language "c++" :nopreprocessoroutput false
            str = Dynare.print(options)
            contains(str, "debug") && contains(str, "nograph") && contains(str, "output=dynamic") && contains(str, "language=c++") && !contains(str, "nopreprocessoroutput")
        catch
            false
        end
    end
end
