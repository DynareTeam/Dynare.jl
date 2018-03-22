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

export @compile, @dynare, @PreprocessorOptions

type PreprocessorOptions
    debug::Bool
    noclearall::Bool
    onlyclearglobals::Bool
    savemacro::Union{Bool,String}
    onlymacro::Bool
    nolinemacro::Bool
    noemptylinemacro::Bool
    notmpterms::Bool
    nolog::Bool
    warn_uninit::Bool
    nograph::Bool
    nointeractive::Bool
    parallel::Union{Bool,String}
    conffile::Union{Bool,String}
    parallel_slave_open_mode::Bool
    parallel_test::Bool
    nostrict::Bool
    stochastic::Bool
    fast::Bool
    minimal_workspace::Bool
    compute_xrefs::Bool
    output::String
    language::String
    params_derivs_order::Int
    json::Union{Bool,String}
    jsonstdout::Bool
    onlyjson::Bool
    jsonderivsimple::Bool
    nopathchange::Bool
    nopreprocessoroutput::Bool
end

"""
    PreprocessorOptions()

Returns default values for the preprocessor options.
"""
function PreprocessorOptions()
    return PreprocessorOptions(false,     # debug
                               false,     # noclearall
                               false,     # onlyclearglobals
                               false,     # savemacro
                               false,     # onlymacro
                               false,     # nolinemacro
                               false,     # noemptylinemacro
                               false,     # notmpterms
                               false,     # nolog
                               false,     # warn_uninit
                               false,     # nograph
                               false,     # nointeractive
                               false,     # parallel
                               false,     # conffile
                               false,     # parallel_slave_open_mode
                               false,     # parallel_test
                               false,     # nostrict
                               false,     # stochastic
                               false,     # fast
                               false,     # minimal_workspace
                               false,     # compute_xrefs
                               "dynamic", # output
                               "julia",   # language
                               0,         # params_derivs_order
                               false,     # json
                               false,     # jsonstdout
                               false,     # onlyjson
                               false,     # jsonderivsimple
                               false,     # nopathchange
                               true)      # nopreprocessoroutput
end


"""
    PreprocessorOptions!(opt::PreprocessorOptions, name::Symbol, value::Any)

Changes the value of option `name` in PreprocessorOptions object `opt`.
"""
function PreprocessorOptions!(opt::PreprocessorOptions, name::Symbol, value::Any)
    if in(name, fieldnames(opt))
        if isa(value, fieldtype(typeof(opt), name))
            setfield!(opt, name, value)
        else
            error("Option `$(string(name))` must be of type $(fieldtype(typeof(opt), name))!")
        end
    else
        error("Option `$(string(name))` is unknown!")
    end
end


"""
    @PreprocessorOptions!(args...)

Sets the preprocessor options (described in the reference manual). If the value of the specified option is not true, the option must be set with a pair (option name and option value).

# Examples
```julia-repl
julia> options = @PreprocessorOptions :nograph :savemacro
Dynare.PreprocessorOptions(false, false, false, true, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, "dynamic", "julia", 0, false, false, false, false, false, true)
```
Sets options `nograph` and `savemacro` to true (default is false).

```julia-repl
julia> options = @PreprocessorOptions :nograph :language "c++" :params_derivs_order 1
Dynare.PreprocessorOptions(false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, "dynamic", "c++", 1, false, false, false, false, false, true)
```
Sets option `nograph` to true, `language` to "c++" (*ie* generates c++ routines instead of julia routines for the model) and `params_derivs_order` to 1 (the preprocessor computes the first order derivates with respect to the parameters).
"""
macro PreprocessorOptions(args...)
    return begin
        n, i = length(eval(args)), 1
        opts = PreprocessorOptions()
        while i<=n
            @assert isa(eval(args[i]), Symbol) "Argument $(string(i)) has to be a Symbol (preprocessor option name)!"
            if i==n || isa(eval(args[i+1]), Symbol)
                PreprocessorOptions!(opts, eval(args[i]), true)
                i+=1
            else
                PreprocessorOptions!(opts, eval(args[i]), eval(args[i+1]))
                i+=2
            end
        end
        opts
    end
end

"""
    print(opt::PreprocessorOptions)

Prints the options values in a string to be passed to the preprocessor binary.

# Examples
```julia-repl
julia> options = @PreprocessorOptions :nograph :savemacro
julia> Dynare.print(options)
" savemacro nograph output=dynamic language=julia nopreprocessoroutput"
```
"""
function print(opt::PreprocessorOptions)
    defaults = PreprocessorOptions()
    optionstring = ""
    for field in fieldnames(defaults)
        if getfield(opt,field)==getfield(defaults,field)
            # User did not change the default value of the preprocessor's option.
            if field==:output
                optionstring = "$optionstring output=dynamic"
            elseif field==:language
                optionstring = "$optionstring language=julia"
            elseif field==:nopreprocessoroutput
                optionstring = "$optionstring nopreprocessoroutput"
            end
        else
            # User changed the value of the preprocessor's option.
            if field==:savemacro
                if typeof(getfield(opt, :savemacro))==String
                    optionstring = "$optionstring savemacro=$(getfield(opt,:savemacro))"
                else # option value is a boolean (true)
                    optionstring = "$optionstring savemacro"
                end
                continue
            end
            if field==:parallel
                if typeof(getfield(opt,:parallel))==String
                    optionstring = "$optionstring parallel=$(getfield(opt,:parallel))"
                else # option value is a boolean (true)
                    optionstring = "$optionstring parallel"
                end
                continue
            end
            if field==:conffile && typeof(getfield(opt,:conffile))==String
                optionstring = "$optionstring conffile=$(getfield(opt,:conffile))"
                continue
            end
            if field==:json
                if typeof(getfield(opt,:json))==String
                    optionstring = "$optionstring json=$(getfield(opt,:json))"
                    continue
                else
                    error("Option `json` must be of type String!")
                end
            end
            if field==:language
                if getfield(opt,:language)!="matlab"
                    optionstring = "$optionstring language=$(getfield(opt,:language))"
                end
                continue
            end
            if field==:output
                optionstring = "$optionstring output=$(getfield(opt,:output))"
                continue
            end
            if field==:params_derivs_order
                if getfield(opt,:params_derivs_order)==1 || getfield(opt,:params_derivs_order)==2
                    optionstring = "$optionstring params_derivs_order=$(getfield(opt,:params_derivs_order))"
                    continue
                else
                    error("Option `params_derivs_order` must be equal to 1 or 2!")
                end
            end
            if getfield(opt, field)==true
                if field!=:nopreprocessoroutput
                    optionstring = "$optionstring $(string(field))"
                end
            end
        end
    end
    return optionstring
end

"""
    compile(modfile::AbstractString, options::PreprocessorOptions)

Compiles `modfile` with `options`, calling the Dynare's preprocessor.
"""
function compile(modfile::AbstractString, options::PreprocessorOptions)
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
    opts = split(print(options))
    run(`$dynare $modfile $opts`)
end

"""
    compile(modfile::AbstractString, options::PreprocessorOptions)

Compiles `modfile` with default options, calling the Dynare's preprocessor.
"""
function compile(modfile::AbstractString)
    opts = PreprocessorOptions()
    compile(modfile, opts)
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

"""
    compile(modfiles::Expr, options::Expr)

Compiles a set of mod files (`modfile` is a vector of strings for the names of the mod files) with default options set with @PreprocessorOptions macro.
"""
macro compile(modfiles::Expr, options::Expr)
    if modfiles.head==:vect && options.head==:macrocall
        for i=1:length(modfiles.args)
            compile(modfiles.args[i], eval(options))
        end
    end
end

"""
    compile(modfiles::AbstractString, options::Expr)

Compiles a mod file (`modfile` is a string for the name of the mod files with or without extension), options are set with the @PreprocessorOptions macro.
"""
macro compile(modfile::AbstractString, options::Expr)
    compile(modfile, eval(options))
end

"""
    compile(modfiles::AbstractString)

Compiles a set of mod files (`modfile` is a vector of strings for the names of the mod files with or without extensions), default options are used.
"""
macro compile(modfiles::Expr)
    options = PreprocessorOptions()
    if modfiles.head==:vect
        for i=1:length(modfiles.args)
            compile(modfiles.args[i], options)
        end
    end
end

"""
    compile(modfiles::AbstractString)

Compiles a mod file (`modfile` is a string for the name of the mod files with or without extension), default options are used.
"""
macro compile(modfile::AbstractString)
    options = PreprocessorOptions()
    compile(modfile, options)
end

function isdynarefile(ext::String)
    if ext==".mod" || ext==".dyn"
        return true
    else
        return false
    end
end

end
