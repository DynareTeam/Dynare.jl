function decision_rules(m::Model, calib::Dict{Symbol, Float64}, steady_state::Vector{Float64})
    @assert length(steady_state) == m.n_endo
    n_dynvars = m.n_back_mixed + m.n_endo + m.n_fwrd_mixed
    yy = [ steady_state[m.zeta_back_mixed]; steady_state; steady_state[m.zeta_fwrd_mixed] ]
    x = zeros(m.n_exo)
    p = calib2vec(m, calib)
    r = Array{Float64}(m.n_endo)
    m.dynamic_mf!(yy, x, p, r)
    jacob = Array{Float64}(m.n_endo, n_dynvars + m.n_exo)
    m.dynamic_mg!(yy, x, p, jacob)

    A = jacob[:, 1:n_dynvars]

    # Deal with static variables
    if m.n_static > 0
        S = jacob[:, m.n_back_mixed + m.zeta_static]
        F = qrfact(S)
        @assert rank(F[:R]) == m.n_static
        A = F[:Q]'*A
    end

    D = [ A[m.n_static+1:end, m.n_back_mixed + m.zeta_back_mixed] A[m.n_static+1:end, m.n_back_mixed+m.n_endo+1:end];
         eye(m.n_back_mixed, m.n_back_mixed)[m.beta_back, :] zeros(m.n_mixed, m.n_fwrd_mixed) ]

    Atilde0plus = A[m.n_static+1:end, m.n_back_mixed + m.zeta_fwrd_mixed]
    Atilde0plus[:, m.beta_fwrd] = zeros(m.n_dynamic, m.n_mixed)
    E = [ -A[m.n_static+1:end, 1:m.n_back_mixed] -Atilde0plus;
         zeros(m.n_mixed, m.n_back_mixed) eye(m.n_fwrd_mixed, m.n_fwrd_mixed)[m.beta_back, :] ]

    (Z, sdim, eigs) = qz!(E, D)

    # BK order condition
    @assert sdim == m.n_back_mixed

    # Compute DR wrt backward endogenous
    Zp = Z'
    Z22 = Zp[m.n_back_mixed+1:end, m.n_back_mixed+1:end]
    Z21 = Zp[m.n_back_mixed+1:end, 1:m.n_back_mixed]

    gyplus = -Z22\Z21

    Z11p = Z[1:m.n_back_mixed, 1:m.n_back_mixed]
    T11 = D[1:m.n_back_mixed, 1:m.n_back_mixed]
    S11 = E[1:m.n_back_mixed, 1:m.n_back_mixed]
    gyminus = Z11p*(T11\S11)/Z11p

    gy = Array{Float64}(m.n_endo, m.n_back_mixed)
    gy[m.zeta_back_mixed, :] = gyminus
    gy[m.zeta_fwrd_mixed, :] = gyplus
    
    if m.n_static > 0
        Abreve0s = A[1:m.n_static, m.n_back_mixed + m.zeta_static]
        Abreve0d = A[1:m.n_static, m.n_back_mixed + m.zeta_dynamic]
        Abreveminus = A[1:m.n_static, 1:m.n_back_mixed]
        Abreveplus = A[1:m.n_static, m.n_back_mixed+m.n_endo+1:end]
        gy[m.zeta_static, :] = -Abreve0s \ (Abreveplus*gyplus*gyminus + Abreve0d*gy[m.zeta_dynamic, :] + Abreveminus)
    end

    # Compute DR wrt shocks
    B = jacob[:, m.n_back_mixed + (1:m.n_endo)]
    B[:, m.zeta_back_mixed] += jacob[:, m.n_back_mixed+m.n_endo+(1:m.n_fwrd_mixed)]*gyplus
    gu = -B\jacob[:,m.n_back_mixed+m.n_endo+m.n_fwrd_mixed+1:end]

    (gy, gu, eigs)
end

function print_decision_rules(m::Model, gy::Matrix{Float64}, gu::Matrix{Float64})
    println("Decision Rules:")
    @printf "%10s " ""
    for j in 1:m.n_endo
        @printf "%10s " m.endo[j]
    end
    @printf "\n"

    for i in 1:m.n_back_mixed
        @printf "%6s(-1) " m.endo[m.zeta_back_mixed[i]]

        for j in 1:m.n_endo
            @printf "%10.6f " gy[j, i]
        end
        @printf "\n"
    end

    for i in 1:m.n_exo
        @printf "%10s " m.exo[i]

        for j in 1:m.n_endo
            @printf "%10.6f " gu[j, i]
        end
        @printf "\n"
    end
end
