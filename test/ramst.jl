# Translation of ramst.mod
#
# Note that we had to introduce an endogenous variable x1 equal to the
# exogenous variable x, because leads on exogenous variables are not yet supported

m =
@modfile begin
    @var c k x1
    @varexo x
    @parameters alph gam delt bet aa
    @model begin
        c + k - aa*x*k(-1)^alph - (1-delt)*k(-1)
        c^(-gam) - (1+bet)^(-1)*(aa*alph*x1(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam)
        x1 = x
    end
end

compute_model_info(m)

calib = Dict(
         :alph => 0.5,
         :gam => 0.5,
         :delt => 0.02,
         :bet => 0.05,
         :aa => 0.5
        )

exoval = Dict( :x => 1.0 )

initval = Dict{Symbol, Float64}()
initval[:k] = ((calib[:delt]+calib[:bet])/(exoval[:x]*calib[:aa]*calib[:alph]))^(1/(calib[:alph]-1))
initval[:c] = calib[:aa]*initval[:k]^calib[:alph]-calib[:delt]*initval[:k]
initval[:x1] = exoval[:x]

# (re)Set initial condition far from the steady state.
initval[:c] = .5*initval[:c]
initval[:k] = .1*initval[:k]

# Solve the static model for the steady state.
s = steady_state(m, calib, initval, exoval)

print_steady_state(m, s)

# Compute a 200 periods perfect foresight simulation
# Start and end at the steady state corresponding to x=1
# There is a shock x=1.2 in t=1
# The results will be in endopath (which includes initial and terminal
#  condition, i.e. it goes from t=0 to t=201)

T = 200
endopath = repmat(s, 1, T+2)
exopath = ones(m.n_exo, T)
exopath[1, 1] = 1.2

perfect_foresight_simul!(m, endopath, exopath, calib)
