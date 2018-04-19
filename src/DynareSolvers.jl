module DynareSolvers

##
 # Copyright (C) 2018 Dynare Team
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

export trustregion

const p1 = .1
const p5 = .5
const p001 = .001
const p0001 = .0001
const macheps = eps(Float64)

type TrustRegionWA
    x::Vector{Float64}                   # Vector of unknowns
    xx::Vector{Float64}                  # Vector of unknowns
    fval0::Vector{Float64}               # residuals of the non linear equations
    fval1::Vector{Float64}               # residuals of the non linear equations (next)
    fjac::Matrix{Float64}                # jacobian matrix of the non linear equations
    fjaccnorm::Vector{Float64}           # norms of the columns of the Jacobian matrix
    fjaccnorm__::Vector{Float64}         # copy of fjaccnorm
    w::Vector{Float64}                   # working array
    p::Vector{Float64}                   # direction
    s::Vector{Float64}
end

function TrustRegionWA(n::Int)
    TrustRegionWA(Vector{Float64}(n), Vector{Float64}(n), Vector{Float64}(n), Vector{Float64}(n), Matrix{Float64}(n,n),
                  Vector{Float64}(n), Vector{Float64}(n), Vector{Float64}(n), Vector{Float64}(n), Vector{Float64}(n))
end


"""
    dogleg!(x::Vector{Float64}, r::Matrix{Float64}, b::Vector{Float64}, d::Vector{Float64}, δ::Float64)

Given an `n` by `n` matrix `r`, an `n` by 1 vector `d` with non zero entries, an `n` by `1`
vector `b`, and a positive number δ, the problem is to determine the convex combination `x`
of the gauss-newton and scaled gradient directions that minimizes (r*x - b) in the least
squares sense, subject to the restriction that the euclidean norm of d*x be at most delta.
"""
function dogleg!(x::Vector{Float64}, r::Matrix{Float64}, b::Vector{Float64}, d::Vector{Float64}, δ::Float64, s::Vector{Float64})
    n = length(x)
    # Compute the Gauss-Newton direction.
    x .= r\b
    # Compute norm of scaled x.
    qnorm = zero(Float64)
    @inbounds for i=1:n
        qnorm += (d[i]*x[i])^2
    end
    qnorm = sqrt(qnorm)
    if qnorm<=δ
        # Gauss-Newton direction is acceptable. There is nothing to do here.
    else
        # Gauss-Newton direction is not acceptable…
        # Compute the scale gradient direction and its norm
        gnorm = zero(Float64)
        @inbounds for i=1:n
            s[i] = zero(Float64)
            @inbounds for j=1:n
                s[i] += r[j,i]*b[j]
            end
            s[i] /= d[i]
            gnorm += s[i]*s[i]
        end
        gnorm = sqrt(gnorm)
        # Compute the norm of the scaled gradient direction.
        # gnorm = norm(s)
        if gnorm>0
            # Normalize and rescale → gradient direction.
            temp0 = zero(Float64)
            @inbounds for i = 1:n
                s[i] /= gnorm*d[i]
            end
            temp0 = zero(Float64)
            @inbounds for i=1:n
                temp1 = zero(Float64)
                @inbounds for j=1:n
                    temp1 += r[i,j]*s[j]
                end
                temp0 += temp1*temp1
            end
            temp0 = sqrt(temp0)
            sgnorm = gnorm/(temp0*temp0)
            if sgnorm<=δ
                # The scaled gradient direction is not acceptable…
                # Compute the point along the dogleg at which the
                # quadratic is minimized.
                bnorm = norm(b)
                temp1 = δ/qnorm
                temp2 = sgnorm/δ
                temp0 = bnorm*bnorm*temp2/(gnorm*qnorm)
                temp0 -= temp1*temp2*temp2-sqrt((temp0-temp1)^2+(one(Float64)-temp1*temp1)*(one(Float64)-temp2*temp2))
                α = temp1*(one(Float64)-temp2*temp2)/temp0
            else
                # The scaled gradient direction is acceptable.
                α = zero(Float64)
            end
        else
            # If the norm of the scaled gradient direction is zero.
            α = δ/qnorm
            sgnorm = .0
        end
        # Form the appropriate  convex combination of the Gauss-Newton direction and the
        # scaled gradient direction.
        temp1 = (one(Float64)-α)*min(sgnorm, δ)
        @inbounds for i = 1:n
            x[i] = α*x[i] + temp1*s[i]
        end
    end
end

"""
    trustregion(f!::Function, j!::Function, y0::Vector{Float64}, tolf::Float64, tolx::Float64, maxiter::Int)

Solves a system of nonlinear equations of several variables using a trust region method.
"""
function trustregion(f!::Function, j!::Function, x0::Vector{Float64}, tolx::Float64, tolf::Float64, maxiter::Int)
    wa = TrustRegionWA(length(x0))
    info = trustregion(f!, j!, x0, tolx, tolf, maxiter, wa)
    return wa.x, info
end

"""
    trustregion(f!::Function, j!::Function, y0::Vector{Float64}, tolf::Float64, tolx::Float64, maxiter::Int, wa::TrustRegionWA)

Solves a system of nonlinear equations of several variables using a trust region method. This version requires a last argument of type `TrustRegionWA`
which holds various working arrays. This version of the solver does not instantiate any array. Results of the solver are available in `wa.x` if `info=1`. 
"""
function trustregion(f!::Function, j!::Function, x0::Vector{Float64}, tolx::Float64, tolf::Float64, maxiter::Int, wa::TrustRegionWA)
    n, iter, info = length(x0), 1, 0
    fnorm, fnorm1 = one(Float64), one(Float64)
    wa.x .= copy(x0)
    # Initial evaluation of the residuals (and compute the norm of the residuals)
    try
        f!(wa.fval0, wa.x)
        fnorm = norm(wa.fval0)
    catch
        error("The system of nonlinear equations cannot be evaluated on the initial guess!")
    end
    # Initial evaluation of the Jacobian
    try
        j!(wa.fjac, wa.x)
    catch
        error("The Jacobian of the system of nonlinear equations cannot be evaluated on the initial guess!")
    end
    # Initialize counters
    ncsucc = zero(Int)
    nslow1 = zero(Int)
    nslow2 = zero(Int)
    # Newton iterations
    while iter<=maxiter && info==0
        # Compute columns norm for the Jacobian matrix.
        @inbounds for i=1:n
            wa.fjaccnorm[i] = zero(Float64)
            @inbounds for j = 1:n
                wa.fjaccnorm[i] += wa.fjac[j,i]*wa.fjac[j,i]
            end
            wa.fjaccnorm[i] = sqrt(wa.fjaccnorm[i])
        end
        if iter==1
            # On the first iteration, calculate the norm of the scaled vector of unknwonws x
            # and initialize the step bound δ. Scaling is done according to the norms of the
            # columns of the initial jacobian.
            @inbounds for i = 1:n
                wa.fjaccnorm__[i] = wa.fjaccnorm[i]
            end
            wa.fjaccnorm__[find(abs.(wa.fjaccnorm__).<macheps)] = one(Float64)
            xnorm = zero(Float64)
            @inbounds for i=1:n
                xnorm += (wa.fjaccnorm__[i]*wa.x[i])^2
            end
            δ = sqrt(xnorm)
            if δ<macheps
                δ = one(Float64)
            end
        else
            wa.fjaccnorm__ .= max.(.1*wa.fjaccnorm__, wa.fjaccnorm)
        end
        # Determine the direction p (with trust region model defined in dogleg routine).
        dogleg!(wa.p, wa.fjac, wa.fval0, wa.fjaccnorm__, δ, wa.s)
        @inbounds for i=1:n
            wa.w[i] = wa.fval0[i]
            @inbounds for j = 1:n
                wa.w[i] -= wa.fjac[i,j]*wa.p[j]
            end
        end
        # Compute the norm of p.
        pnorm = zero(Float64)
        @inbounds for i=1:n
            pnorm += (wa.fjaccnorm__[i]*wa.p[i])^2
        end
        pnorm = sqrt(pnorm)
        # On first iteration adjust the initial step bound.
        if iter==1
            δ = min(δ, pnorm)
        end
        # Evaluate the function at xx = x+p and calculate its norm.
        @inbounds for i = 1:n
            wa.xx[i] = wa.x[i]-wa.p[i]
        end
        f!(wa.fval1, wa.xx)
        fnorm1 = norm(wa.fval1)
        # Compute the scaled actual reduction.
        if fnorm1<fnorm
            actualreduction = one(Float64)-(fnorm1/fnorm)^2
        else
            actualreduction = -one(Float64)
        end
        # Compute the scaled predicted reduction and the ratio of the actual to the
        # predicted reduction.
        t = norm(wa.w)
        ratio = zero(Float64)
        if t<fnorm
            predictedreduction = one(Float64) - (t/fnorm)^2
            ratio = actualreduction/predictedreduction
        end
        # Update the radius of the trust region if need be.
        if ratio<p1
            # Reduction is much smaller than predicted… Reduce the radius of the trust region.
            ncsucc = 0
            δ = p5*δ
        else
            ncsucc += 1
            if ratio>=p5 || ncsucc>1
                δ = max(δ,pnorm/p5)
            end
            if abs(ratio-one(Float64))<p1
                δ = pnorm/p5
            end
        end
        if ratio>=1.0e-4
            # Succesfull iteration. Update x, xnorm, fval0, fnorm and fjac.
            xnorm = zero(Float64)
            @inbounds for i=1:n
                wa.x[i] = wa.xx[i]
                xnorm += (wa.fjaccnorm__[i]*wa.x[i])^2
                wa.fval0[i] = wa.fval1[i]
            end
            xnorm = sqrt(xnorm)
            fnorm = fnorm1
        end
        # Determine the progress of the iteration.
        nslow1 += 1
        if actualreduction>=p001
            nslow1 = 0
        end
        # Test for convergence.
        if δ<tolx*xnorm || fnorm<tolf
            info = 1
            continue
        end
        # Tests for termination and stringent tolerances.
        if p1*max(p1*δ, pnorm)<=macheps*xnorm
            # xtol is too small. no further improvement in
            # the approximate solution x is possible.
            info = 3
            continue
        end
        if nslow1==15
            # iteration is not making good progress, as
            # measured by the improvement from the last
            # fifteen iterations.
            info = 4
            continue
        end
        # Compute the jacobian for next the iteration.
        j!(wa.fjac, wa.x)
        iter += 1
    end
    if info==0 && iter>maxiter
        info = 2
        fill!(wa.x, Inf)
    end
    return info
end

end
