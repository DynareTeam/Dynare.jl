# Translation of example1.mod

m =
@modfile begin
    @var y c k a h b
    @varexo e u
    @parameters beta rho alpha delta theta psi tau
    @model begin
        c*theta*h^(1+psi) = (1-alpha)*y
        k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))*(exp(b(+1))*alpha*y(+1)+(1-delta)*k))
        y = exp(a)*(k(-1)^alpha)*(h^(1-alpha))
        k = exp(b)*(y-c)+(1-delta)*k(-1)
        a = rho*a(-1)+tau*b(-1) + e
        b = tau*a(-1)+rho*b(-1) + u
    end
end

compute_model_info(m)

calib = Dict(
         :alpha => 0.36,
         :rho   => 0.95,
         :tau   => 0.025,
         :beta  => 0.99,
         :delta => 0.025,
         :psi   => 0.0,
         :theta => 2.95
        )

initval = Dict(
           :y => 1.08068253095672,
           :c => 0.80359242014163,
           :h => 0.29175631001732,
           :k => 11.08360443260358,
           :a => 0.0,
           :b => 0.0
           )

s = steady_state(m, calib, initval)

print_steady_state(m, s)

(gy, gu, eigs) = decision_rules(m, calib, s)

println("Eigenvalues: ", eigs)
println()

print_decision_rules(m, gy, gu)
