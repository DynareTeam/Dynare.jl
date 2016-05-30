module Dynare

using NLsolve

export @var, @varexo, @parameters, @model, @modfile

typealias MathExpr Union{Expr,Symbol,Number}

type Model
    endo::Vector{Symbol}
    exo::Vector{Symbol}
    param::Vector{Symbol}
    eqs::Vector{MathExpr}

    n_endo::Int
    n_exo::Int
    n_param::Int
    zeta_static::Vector{Int}
    zeta_back::Vector{Int}
    zeta_fwrd::Vector{Int}
    zeta_mixed::Vector{Int}
    zeta_back_mixed::Vector{Int}
    zeta_fwrd_mixed::Vector{Int}
    zeta_dynamic::Vector{Int}
    beta_fwrd::Vector{Int}
    pi_fwrd::Vector{Int}
    beta_back::Vector{Int}
    pi_back::Vector{Int}
    n_static::Int
    n_back::Int
    n_fwrd::Int
    n_mixed::Int
    n_back_mixed::Int
    n_fwrd_mixed::Int
    n_dynamic::Int
    static_mf!::Function # Static model function (y,x,p,residual)
    static_mg!::Function # Static model jacobian (y,x,p,jacobian)
    dynamic_mf!::Function # Dynamic model function (yy,x,p,residual) where yy = (y-, y, y+)
    dynamic_mg!::Function # Dynamic model Jacobian (yy,x,p,jacobian)
end

Model() = Model(Array(Symbol, 0), Array(Symbol, 0), Array(Symbol, 0), Array(MathExpr, 0),
                0, 0, 0,
                Array(Int, 0), Array(Int, 0),
                Array(Int, 0), Array(Int, 0),
                Array(Int, 0), Array(Int, 0),
                Array(Int, 0),
                Array(Int, 0), Array(Int, 0),
                Array(Int, 0), Array(Int, 0),
                0, 0, 0, 0, 0, 0, 0,
                (y::Vector{Float64}, x::Vector{Float64}, p::Vector{Float64}, fvec::Vector{Float64}) -> nothing,
                (y::Vector{Float64}, x::Vector{Float64}, p::Vector{Float64}, fjac::Matrix{Float64}) -> nothing,
                (yy::Vector{Float64}, x::Vector{Float64}, p::Vector{Float64}, fvec::Vector{Float64}) -> nothing,
                (yy::Vector{Float64}, x::Vector{Float64}, p::Vector{Float64}, fjac::Matrix{Float64}) -> nothing)

macro var(x...)
    M_.endo = [x...]
end


macro varexo(x...)
    M_.exo = [x...]
end

macro parameters(x...)
    M_.param = [x...]
end

# Function to remove some :line expressions that appear sometimes
remove_line_blocks(e::Any) = e

function remove_line_blocks(e::Expr)
    if e.head == :block
        assert(sum(map(y->(typeof(y) == Symbol || y.head != :line), e.args)) == 1)
        for x in e.args
            if typeof(x) == Symbol || x.head != :line
                return x
            end
        end
    else
        e.args = map(remove_line_blocks, e.args)
        return e
    end
end

macro model(x)
    M_.eqs = map(remove_line_blocks, filter(y->(typeof(y)!=LineNumberNode && (typeof(y) != Expr || y.head != :line)), x.args))
end

macro modfile(x)
    global M_ = Model()
    eval(x)
    M_
end

# Fills a dictionnary which associates to each symbol its max lead and max lag
scan_expr(e::Number, endo::Vector{Symbol}, lead_lags::Dict{Symbol,Tuple{Int,Int}}) = nothing
scan_expr(e::Symbol, endo::Vector{Symbol}, lead_lags::Dict{Symbol,Tuple{Int,Int}}) = nothing

function scan_expr(e::Expr, endo::Vector{Symbol}, lag_lead::Dict{Symbol,Tuple{Int,Int}})
    if e.head == :call
        if !isempty(filter(x->x==e.args[1], endo))
            v = e.args[1]
            shift = e.args[2]
            @assert isa(shift, Integer)
            (la, le) = lag_lead[v]
            la = max(-shift, la)
            le = max(shift, le)
            lag_lead[v] = (la, le)
        else
            map(x->scan_expr(x,endo,lag_lead), e.args[2:end])
        end
    elseif e.head == :(=)
        map(x->scan_expr(x,endo,lag_lead), e.args)
    else
        error("Unknown head")
    end
end

# Computes the zeta_* and their lengths (the n_*)
function compute_var_categories(m::Model)
    d = Dict{Symbol,Tuple{Int,Int}}()
    for i = 1:m.n_endo
        d[m.endo[i]] = (0,0)
    end
    for j = 1:m.n_endo
        scan_expr(m.eqs[j], m.endo, d)
    end
    for i = 1:m.n_endo
        (la,le) = d[m.endo[i]]
        if la > 1 || le > 1
            error("Lags or leads of more than 1 are not supported")
        elseif la == 0 && le == 0
            push!(m.zeta_static, i)
        elseif la == 0 && le == 1
            push!(m.zeta_fwrd, i)
        elseif la == 1 && le == 0
            push!(m.zeta_back, i)
        else
            push!(m.zeta_mixed, i)
        end
    end
    m.zeta_back_mixed = sort!(union(m.zeta_back, m.zeta_mixed))
    m.zeta_fwrd_mixed = sort!(union(m.zeta_fwrd, m.zeta_mixed))
    m.zeta_dynamic = setdiff(collect(1:m.n_endo), m.zeta_static)

    # Compute the length of the 7 vectors
    for t = (:back, :fwrd, :static, :mixed, :back_mixed, :fwrd_mixed, :dynamic)
        nt = symbol(string("n_", t))
        zt = symbol(string("zeta_", t))
        setfield!(m, nt, length(getfield(m, zt)))
    end

    # Compute beta_back and beta_fwrd
    for i = 1:m.n_fwrd_mixed
        if m.zeta_fwrd_mixed[i] in m.zeta_mixed
            push!(m.beta_fwrd, i)
        else
            push!(m.pi_fwrd, i)
        end
    end

    for i = 1:m.n_back_mixed
        if m.zeta_back_mixed[i] in m.zeta_mixed
            push!(m.beta_back, i)
        else
            push!(m.pi_back, i)
        end
    end
end

import Calculus.differentiate

include("symbolic.jl")

function compute_static_mf(m::Model)
    reqs = Array(MathExpr, m.n_endo)
    jacob = Array(MathExpr, m.n_endo, m.n_endo)

    # Compute expressions for static residuals and jacobian
    for i = 1:m.n_endo
        reqs[i] = to_static(homogeneize(m.eqs[i]), m.endo)
        for j = 1:m.n_endo
            jacob[i,j] = differentiate(reqs[i], m.endo[j])
        end
    end

    # Prepare the exprs for incorporation in static_mf
    env = Dict{Symbol,MathExpr}()
    for i = 1:m.n_endo
        env[m.endo[i]] = :(y[$i])
    end
    for i = 1:m.n_exo
        env[m.exo[i]] = :(x[$i])
    end
    for i = 1:m.n_param
        env[m.param[i]] = :(p[$i])
    end
    for i = 1:m.n_endo
        reqs[i] = subst_symb(reqs[i], env)
        for j = 1:m.n_endo
            jacob[i,j] = subst_symb(jacob[i,j], env)
        end
    end

    # Construct static model functions
    f_assigns = Array(Expr, 0)
    g_assigns = Array(Expr, 0)
    for i = 1:m.n_endo
        push!(f_assigns, :(fvec[$i] = $(reqs[i])))
        for j = 1:m.n_endo
            push!(g_assigns, :(fjac[$i, $j] = $(jacob[i,j])))
        end
    end

    @eval begin
        function static_mf!(y::Vector{Float64}, x::Vector{Float64}, p::Vector{Float64}, fvec::Vector{Float64})
            $(f_assigns...)
        end
    end

    @eval begin
        function static_mg!(y::Vector{Float64}, x::Vector{Float64}, p::Vector{Float64}, fjac::Matrix{Float64})
            $(g_assigns...)
        end
    end

    m.static_mf! = static_mf!
    m.static_mg! = static_mg!
end

function compute_dynamic_mf(m::Model)
    ndyn = m.n_back_mixed + m.n_endo + m.n_fwrd_mixed
    njcols = ndyn + m.n_exo
    reqs = Array(MathExpr, m.n_endo)
    jacob = Array(MathExpr, m.n_endo, njcols)

    # Create gensyms for lagged and leaded variables
    lag_gensyms = Array(Symbol, m.n_back_mixed)
    lag_dict = Dict{Symbol,Symbol}()
    for i = 1:m.n_back_mixed
        lag_gensyms[i] = gensym()
        lag_dict[m.endo[m.zeta_back_mixed[i]]] = lag_gensyms[i]
    end
    lead_gensyms = Array(Symbol, m.n_fwrd_mixed)
    lead_dict = Dict{Symbol,Symbol}()
    for i = 1:m.n_fwrd_mixed
        lead_gensyms[i] = gensym()
        lead_dict[m.endo[m.zeta_fwrd_mixed[i]]] = lead_gensyms[i]
    end

    dynvars = [ lag_gensyms; m.endo; lead_gensyms ]
    
    # Compute expressions for static residuals and jacobian
    for i = 1:m.n_endo
        reqs[i] = subst_lag_lead(homogeneize(m.eqs[i]), lag_dict, lead_dict)
        for j = 1:ndyn
            jacob[i,j] = differentiate(reqs[i], dynvars[j])
        end
        for j = 1:m.n_exo
            jacob[i,ndyn+j] = differentiate(reqs[i], m.exo[j])
        end
    end

    # Prepare the exprs for incorporation in static_mf
    env = Dict{Symbol,MathExpr}()
    for i = 1:ndyn
        env[dynvars[i]] = :(yy[$i])
    end
    for i = 1:m.n_exo
        env[m.exo[i]] = :(x[$i])
    end
    for i = 1:m.n_param
        env[m.param[i]] = :(p[$i])
    end
    for i = 1:m.n_endo
        reqs[i] = subst_symb(reqs[i], env)
        for j = 1:njcols
            jacob[i,j] = subst_symb(jacob[i,j], env)
        end
    end
    
    # Construct dynamic model functions
    f_assigns = Array(Expr, 0)
    g_assigns = Array(Expr, 0)
    for i = 1:m.n_endo
        push!(f_assigns, :(fvec[$i] = $(reqs[i])))
        for j = 1:njcols
            push!(g_assigns, :(fjac[$i, $j] = $(jacob[i,j])))
        end
    end

    @eval begin
        function dynamic_mf!(yy::Vector{Float64}, x::Vector{Float64}, p::Vector{Float64}, fvec::Vector{Float64})
            $(f_assigns...)
        end
    end

    @eval begin
        function dynamic_mg!(yy::Vector{Float64}, x::Vector{Float64}, p::Vector{Float64}, fjac::Matrix{Float64})
            $(g_assigns...)
        end
    end

    m.dynamic_mf! = dynamic_mf!
    m.dynamic_mg! = dynamic_mg!
end

export compute_model_info

function compute_model_info(m::Model)
    m.n_endo = length(m.endo)
    @assert m.n_endo == length(m.eqs)
    m.n_exo = length(m.exo)
    m.n_param = length(m.param)
    compute_var_categories(m)
    compute_static_mf(m)
    compute_dynamic_mf(m)
    nothing
end

function calib2vec(m::Model, calib::Dict{Symbol, Float64})
    p = Array(Float64, m.n_param)
    for i in 1:m.n_param
        p[i] = calib[m.param[i]]
    end
    return(p)
end

function initval2vec(m::Model, initval::Dict{Symbol, Float64})
    y = Array(Float64, m.n_endo)
    for i in 1:m.n_endo
        y[i] = initval[m.endo[i]]
    end
    return(y)
end

function exoval2vec(m::Model, exoval::Dict{Symbol, Float64})
    x = Array(Float64, m.n_exo)
    for i in 1:m.n_exo
        x[i] = exoval[m.exo[i]]
    end
    return(x)
end

function steady_state(m::Model, calib::Dict{Symbol, Float64}, initval::Dict{Symbol, Float64}, exoval::Dict{Symbol, Float64})
    p = calib2vec(m, calib)
    iv = initval2vec(m, initval)
    ev = exoval2vec(m, exoval)
    mf!(y::Vector{Float64}, fvec::Vector{Float64}) = m.static_mf!(y, ev, p, fvec)
    mg!(y::Vector{Float64}, fjac::Matrix{Float64}) = m.static_mg!(y, ev, p, fjac)
    r = nlsolve(mf!, mg!, iv, show_trace = true)
    @assert converged(r)
    return(r.zero)
end

# Version without exogenous variables, which are assumed to have a zero value
function steady_state(m::Model, calib::Dict{Symbol, Float64}, initval::Dict{Symbol, Float64})
    exoval = Dict{Symbol, Float64}()
    for i in 1:m.n_exo
        exoval[m.exo[i]] = 0.0
    end
    return steady_state(m, calib, initval, exoval)
end

function print_steady_state(m::Model, steady_state::Vector{Float64})
    println("Steady State:")
    for i in 1:m.n_endo
        @printf "%10s %10.6f\n" m.endo[i] steady_state[i]
    end
    println()
end

export steady_state, print_steady_state

include("qz.jl")
include("decision_rules.jl")

export decision_rules, print_decision_rules

include("perfect_foresight_simul.jl")

export perfect_foresight_simul!

end
