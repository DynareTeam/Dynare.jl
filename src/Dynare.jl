module Dynare

##
 # Copyright (C) 2015-2018 Dynare Team
 #
 # This file is part of Dynare.
 #
 # Dynare is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # Dynare is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
##

# Check that the preprocessor is here where it should be...
const depsjlpath = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")

if !isfile(depsjlpath)
    error("Dynare's preprocessor not installed properly, run Pkg.build(\"Dynare\"), restart Julia and try again.")
else
    include(depsjlpath)
end

export @compile, @dynare

function compile(modfile::AbstractString)
    # Add current path to LOAD_PATH if necessary.
    if isempty(findin([pwd()], LOAD_PATH))
        unshift!(LOAD_PATH, pwd())
    end
    # Append extension if necessary and check extension.
    basename, ext = splitext(modfile)
    if isempty(ext)
        modfile = "$modfile.mod"
    else
        if ~isdynarefile(ext)
            error("The Dynare model file must have a mod or dyn extension!")
        end
    end
    # Call the preprocessor.
    run(`$dynare $modfile language=julia output=dynamic nopreprocessoroutput`)
end

macro dynare(modfiles...)
    ex = Expr(:toplevel)
    if length(modfiles)>1
        for modfile in modfiles
            eval(:(compile($modfile)))
            basename = split(modfile, ".mod"; keep=false)
            push!(ex.args, Expr(:import, Symbol(basename[1])))
        end
    else
        eval(:(compile($modfiles)))
        basename = split(modfiles[1], ".mod"; keep=false)
        push!(ex.args, Expr(:importall, Symbol(basename[1])))
    end
    return ex
end

macro compile(modfiles::Expr)
    if modfiles.head==:vect
        for i=1:length(modfiles.args)
            compile(modfiles.args[i])
        end
    end
end

macro compile(modfile::AbstractString)
    #return :(compile($modfile))
    compile(modfile)
end

function isdynarefile(ext::String)
    if ext==".mod" || ext==".dyn"
        return true
    else
        return false
    end
end

end
