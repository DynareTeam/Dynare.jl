module TestFunctions

import Base.sign

export rosenbrock, rosenbrock!, powell1, powell1!, powell2, powell2!, wood, wood!
export helicalvalley, helicalvalley!, watson, watson!, chebyquad, chebyquad!
export brown, brown!, discreteboundaryvalue, discreteboundaryvalue!
export discreteintegralequation, discreteintegralequation!
export trigonometric, trigonometric!, variablydimensioned, variablydimensioned!
export broydentridiagonal, broydentridiagonal!
export broydenbanded, broydenbanded! 

const ZERO = 0.0
const ONE = 1.0
const TWO = 2.0
const THREE = 3.0
const FOUR = 4.0
const FIVE = 5.0
const SIX = 6.0
const EIGHT = 8.0
const TEN = 10.0
const FIFTN = 15.0
const TWENTY = 20.0
const HUNDRD = 100.0
const HALF = 0.5
const C1 = 1.0e4
const C2 = 1.0001
const C3 = 2.0e2
const C4 = 20.2
const C5 = 19.8
const C6 = 180.0
const C7 = 0.25
const C8 = 0.5
const C9 = 29.0
const TPI = EIGHT*atan(ONE)

"""
    sign(a::Real, b::Real)

Returns |a|*sign(b).
"""
function sign(a::Real, b::Real)
    return abs(a)*sign(b)
end

"""
    rosenbrock!(fval::Vector{Float64}, x::Vector{Float64})

Evaluates the Rosenbrock function. 
"""
function rosenbrock!(fval::Vector{Float64}, x::Vector{Float64})
    fval[1] = ONE-x[1]
    fval[2] = TEN*(x[2]-fval[1]*fval[1])
end

"""
    rosenbrock!(fval::Matrix{Float64}, x::Vector{Float64})

Evaluates the jacobian matrix of the Rosenbrock function. 
"""
function rosenbrock!(fjac::Matrix{Float64}, x::Vector{Float64})
    fjac[1,1] = -ONE
    fjac[1,2] = ZERO
    fjac[2,1] = -TWENTY*x[1]
    fjac[2,2] = TEN
end

"""
    rosenbrock()

Provides initial guess for solving the nonlinear system of equation based on the Rosenbrock function. 
"""
function rosenbrock()
    return [-C1, ONE]
end

function powell1!(fval::Vector{Float64}, x::Vector{Float64})
    fval[1] = x[1]+TEN*x[2]
    fval[2] = sqrt(FIVE)*(x[3]-x[4])
    tmp = x[2]-TWO*x[3]
    fval[3] = tmp*tmp
    tmp = x[1]-x[4]
    fval[4] = sqrt(TEN)*tmp*tmp
end

function powell1!(fjac::Matrix{Float64}, x::Vector{Float64})
    fill!(fjac, 0.0)
    fjac[1,1] = ONE
    fjac[1,2] = TEN
    fjac[2,3] = sqrt(FIVE)
    fjac[2,4] = -fjac[2,3]
    fjac[3,2] = TWO*(x[2]-TWO*x[3])
    fjac[3,3] = -TWO*fjac[3,2]
    fjac[4,1] = TWO*sqrt(TEN)*(x[1]-x[4])
    fjac[4,4] = -fjac[4,1]
end

function powell1()
    return [THREE, -ONE, ZERO, ONE]
end

function powell2!(fval::Vector{Float64}, x::Vector{Float64})
    fval[1] = C1*x[1]*x[2]-ONE
    fval[2] = exp(-x[1])+exp(-x[2])-C2
end

function powell2!(fjac::Matrix{Float64}, x::Vector{Float64})
    fjac[1,1] = C1*x[2]
    fjac[1,2] = C1*x[1]
    fjac[2,1] = -exp(-x[1])
    fjac[2,2] = -exp(-x[2])
end

function powell2()
    return [ZERO, ONE]
end

function wood!(fval::Vector{Float64}, x::Vector{Float64})
    tmp1 = x[2]-x[1]*x[1]
    tmp2 = x[4]-x[3]*x[3]
    fval[1] = -C3*x[1]*tmp1-(ONE-x[1])
    fval[2] = C3*tmp1+C4*(x[2]-1.0)+C5*(x[4]-ONE)
    fval[3] = -C6*x[3]*tmp2-(ONE-x[3])
    fval[4] = C6*tmp2+C4*(x[4]-ONE)+C5*(x[2]-ONE)
end

function wood!(fjac::Matrix{Float64}, x::Vector{Float64})
    fill!(fjac, ZERO)
    tmp1 = x[2]-THREE*x[1]*x[1]
    tmp2 = x[4]-THREE*x[3]*x[3]
    fjac[1,1] = -C3*tmp1+ONE
    fjac[1,2] = -C3*x[1]
    fjac[2,1] = -TWO*C3*x[1]
    fjac[2,2] = C3+C4
    fjac[2,4] = C5
    fjac[3,3] = -C6*tmp2+ONE
    fjac[3,4] = -C6*x[3]
    fjac[4,2] = C5
    fjac[4,3] = -TWO*C6*x[3]
    fjac[4,4] = C6+C4
end

function wood()
    return [-THREE, -ONE, -THREE, -ONE]
end

function helicalvalley!(fval::Vector{Float64}, x::Vector{Float64})
    tmp1 = sign(C7, x[2])
    if x[1]>ZERO
        tmp1 = atan(x[2]/x[1])/TPI
    end
    if x[1]<ZERO
        tmp1 = atan(x[2]/x[1])/TPI+C8
    end
    tmp2 = sqrt(x[1]*x[1]+x[2]*x[2])
    fval[1] = TEN*(x[3]-TEN*tmp1)
    fval[2] = TEN*(tmp2-ONE)
    fval[3] = x[3]
end

function helicalvalley!(fjac::Matrix{Float64}, x::Vector{Float64})
    tmp = x[1]*x[1]+x[2]*x[2]
    tmp1 = TPI*tmp
    tmp2 = sqrt(tmp)
    fjac[1,1] = HUNDRD*x[2]/tmp1
    fjac[1,2] = -HUNDRD*x[1]/tmp1
    fjac[1,3] = TEN
    fjac[2,1] = TEN*x[1]/tmp2
    fjac[2,2] = TEN*x[2]/tmp2
    fjac[2,3] = ZERO
    fjac[3,1] = ZERO
    fjac[3,2] = ZERO
    fjac[3,3] = ONE
end

function helicalvalley()
    return [-ONE, ZERO, ZERO]
end

function watson!(fval::Vector{Float64}, x::Vector{Float64})
    n = length(fval)
    fill!(fval, ZERO)
    for i=1:29
        ti = Float64(i)/C9
        sum1 = ZERO
        tmp = ONE
        for j=2:n
            sum1 += Float64(j-1)*tmp*x[j]
            tmp *= ti
        end
        sum2 = ZERO
        tmp = ONE
        for j=1:n
            sum2 += tmp*x[j]
            tmp *= ti
        end
        tmp1 = sum1-sum2*sum2-ONE
        tmp2 = TWO*ti*sum2
        tmp = ONE/ti
        for k=1:n
            fval[k] += tmp*(Float64(k-1)-tmp2)*tmp1
            tmp *= ti
        end
    end
    tmp = x[2]-x[1]*x[1]-ONE
    fval[1] += x[1]*(ONE-TWO*tmp)
    fval[2] += tmp
end

function watson!(fjac::Matrix{Float64}, x::Vector{Float64})
    n = length(x)
    fill!(fjac, ZERO)
    for i=1:29
        ti = Float64(i)/C9
        sum1 = ZERO
        tmp = ONE
        for j=2:n
            sum1 += Float64(j-1)*tmp*x[j]
            tmp *= ti 
        end
        sum2 = ZERO
        tmp = ONE
        for j=1:n
            sum2 += tmp*x[j]
            tmp *= ti
        end
        tmp1 = TWO*(sum1-sum2*sum2-ONE)
        tmp2 = TWO*sum2
        tmp = ti*ti
        tk = ONE
        for k=1:n
            tj = tk
            for j=k:n
                fjac[k,j] += tj*((Float64(k-1)/ti-tmp2)*(Float64(j-1)/ti-tmp2)-tmp1)
                tj *= ti
            end
            tk *= tmp
        end
    end
    fjac[1,1] = fjac[1,1]+SIX*x[1]*x[1]-TWO*x[2]+THREE
    fjac[1,2] = fjac[1,2]-TWO*x[1]
    fjac[2,2] = fjac[2,2]+ONE
    for k=1:n
        for j=k:n
            fjac[j,k] = fjac[k,j]
        end
    end
end

function watson(n::Int)
    return zeros(Float64, n)
end

function chebyquad!(fval::Vector{Float64}, x::Vector{Float64})
    n = length(fval)
    fill!(fval, ZERO)
    for j=1:n
        tmp1 = ONE
        tmp2 = TWO*x[j]-ONE
        tmp = TWO*tmp2
        for i=1:n
            fval[i] += tmp2
            ti = tmp*tmp2-tmp1
            tmp1 = tmp2
            tmp2 = ti
        end
    end
    tk = ONE/Float64(n)
    iev = -1
    for k=1:n
        fval[k] *= tk
        if iev>0
            fval[k] += ONE/(Float64(k)^2-ONE)
        end
        iev = -iev
    end
end

function chebyquad!(fjac::Matrix{Float64}, x::Vector{Float64})
    n = length(x)
    tk = ONE/Float64(n)
    for j=1:n
        tmp1 = ONE
        tmp2 = TWO*x[j]-ONE
        tmp = TWO*tmp2
        tmp3 = ZERO
        tmp4 = TWO
        for k=1:n
            fjac[k,j] = tk*tmp4
            ti = FOUR*tmp2+tmp*tmp4-tmp3
            tmp3 = tmp4
            tmp4 = ti
            ti = tmp*tmp2-tmp1
            tmp1 = tmp2
            Tmp2 = ti
        end
    end
end

function chebyquad(n::Int)
    h = ONE/Float64(n+1)
    x = Vector{Float64}(n)
    for j=1:n
        x[j] = Float64(j)*h
    end
    return x
end

function brown!(fval::Vector{Float64}, x::Vector{Float64})
    n = length(fval)
    sum = -Float64(n+1)
    prod = ONE
    for j=1:n
        sum += x[j]
        prod *= x[j]
    end
    for k=1:n-1
        fval[k]=x[k]+sum
    end
    fval[n] = prod-ONE
end

function brown!(fjac::Matrix{Float64}, x::Vector{Float64})
    n = length(x)
    prod = ONE
    for j=1:n
        prod *= x[j]
        for k=1:n
            fjac[k,j] = ONE
        end
        fjac[j,j] = TWO
    end
    for j=1:n
        tmp = x[j]
        if abs(tmp)<eps(Float64)
            tmp = ONE
            prod = ONE
            for k=1:n
                if k!=j
                    prod *= x[k]
                end
            end
        end
        fjac[n,j] = prod/tmp
    end
end

function brown(n::Int)
    x = Vector{Float64}(n)
    fill!(x, HALF)
    return x
end

function discreteboundaryvalue!(fval::Vector{Float64}, x::Vector{Float64})
    n = length(fval)
    h = ONE/Float64(n+1)
    hh = h*h
    for k=1:n
        tmp = (x[k]+Float64(k)*h+ONE)^3
        tmp1 = ZERO
        if k!=1
            tmp1 = x[k-1]
        end
        tmp2 = ZERO
        if k!=n
            tmp2 = x[k+1]
        end
        fval[k] = TWO*x[k]-tmp1-tmp2+tmp*hh/TWO        
    end
end

function discreteboundaryvalue!(fjac::Matrix{Float64}, x::Vector{Float64})
    n = length(x)
    h = ONE/Float64(n+1)
    hh = h*h
    fill!(fjac, ZERO)
    for k=1:n
        tmp = THREE*(x[k]+Float64(k)*h+ONE)^2
        fjac[k,k] = TWO+tmp*hh/TWO
        if k!=1
            fjac[k,k-1] = -ONE
        end
        if k!=n
            fjac[k,k+1] = -ONE
        end
    end
end

function discreteboundaryvalue(n::Int)
    h = ONE/Float64(n+1)
    x = Vector{Float64}(n)
    for j=1:n
        tj = Float64(j)*h
        x[j] = tj*(tj-ONE)
    end
    return x
end

function discreteintegralequation!(fval::Vector{Float64}, x::Vector{Float64})
    n = length(fval)
    h = ONE/Float64(n+1)
    for k=1:n
        tk = Float64(k)*h
        sum1 = ZERO
        for j=1:k
            tj = Float64(j)*h
            tmp = (x[j]+tj+ONE)^3
            sum1 += tj*tmp
        end
        sum2 = ZERO
        kp1 = k+1
        if n>=kp1
            for j=kp1:n
                tj = Float64(j)*h
                tmp = (x[j]+tj+ONE)^3
                sum2 += (ONE-tj)*tmp
            end
        end
        fval[k] = x[k] + h*((ONE-tk)*sum1+tk*sum2)/TWO
    end
end

function discreteintegralequation!(fjac::Matrix{Float64}, x::Vector{Float64})
    n = length(x)
    h = ONE/Float64(n+1)
    for k=1:n
        tk = Float64(k)*h
        for j=1:n
            tj = Float64(j)*h
            tmp = THREE*(x[j]+tj+ONE)^2
            fjac[k,j] = h*min(tj*(ONE-tk),tk*(ONE-tj))*tmp/TWO
        end
        fjac[k,k] += ONE
    end
end

discreteintegralequation(n::Int) = discreteboundaryvalue(n)

function trigonometric!(fval::Vector{Float64}, x::Vector{Float64})
    n = length(fval)
    sum = ZERO
    for j=1:n
        fval[j] = cos(x[j])
        sum += fval[j] 
    end
    for k=1:n
        fval[k] = Float64(n+k)-sin(x[k])-sum-Float64(k)*fval[k]
    end
end

function trigonometric!(fjac::Matrix{Float64}, x::Vector{Float64})
    n = length(x)
    for j=1:n
        tmp = sin(x[j])
        for k=1:n
            fjac[k,j] = tmp
        end
        fjac[j,j] = Float64(j+1)*tmp-cos(x[j])
    end
end

function trigonometric(n::Int)
    x = Vector{Float64}(n)
    fill!(x, ONE/Float64(n))
    return x
end

function variablydimensioned!(fval::Vector{Float64}, x::Vector{Float64})
    n = length(fval)
    sum = ZERO
    for j=1:n
        sum += Float64(j)*(x[j]-ONE)
    end
    tmp = sum*(ONE+TWO*sum*sum)
    for k=1:n
        fval[k] = x[k]-ONE+Float64(k)*tmp
    end
end

function variablydimensioned!(fjac::Matrix{Float64}, x::Vector{Float64})
    n = length(x)
    sum = ZERO
    for j=1:n
        sum += Float64(j)*(x[j]-ONE)
    end
    tmp = ONE+SIX*sum*sum
    for k=1:n
        for j=1:n
            fjac[k,j] = Float64(k*j)*tmp
            fjac[j,k] = fjac[k,j]
        end
        fjac[k,k] += ONE
    end
end

function variablydimensioned(n::Int)
    h = ONE/Float64(n)
    x = Vector{Float64}(n)
    for j=1:n
        x[j] = ONE-Float64(j)*h
    end
    return x
end

function broydentridiagonal!(fval::Vector{Float64}, x::Vector{Float64})
    n = length(fval)
    for k=1:n
        tmp = (THREE-TWO*x[k])*x[k]
        tmp1 = ZERO
        if k!=1
            tmp1 = x[k-1]
        end
        tmp2 = ZERO
        if k!=n
            tmp2 = x[k+1]
        end
        fval[k] = tmp-tmp1-TWO*tmp2+ONE
    end
end

function broydentridiagonal!(fjac::Matrix{Float64}, x::Vector{Float64})
    n = length(x)
    fill!(fjac, ZERO)
    for k=1:n
        fjac[k,k] = THREE-FOUR*x[k]
        if k!=1
            fjac[k,k-1] = -ONE
        end
        if k!=n
            fjac[k,k+1] = -TWO
        end
    end
end

function broydentridiagonal(n::Int)
    x = Vector{Float64}(n)
    fill!(x, -ONE)
    return x
end

function broydenbanded!(fval::Vector{Float64}, x::Vector{Float64})
    n = length(fval)
    ml = 5
    mu = 1
    for k=1:n
        k1 = max(1,k-ml)
        k2 = min(k+mu,n)
        tmp = ZERO
        for j=k1:k2
            if j!=k
                tmp += x[j]*(ONE+x[j])
            end
        end
        fval[k] = x[k]*(TWO+FIVE*x[k]*x[k])+ONE-tmp
    end
end

function broydenbanded!(fjac::Matrix{Float64}, x::Vector{Float64})
    n = length(x)
    ml = 5
    mu = 1
    fill!(fjac,ZERO)
    for k=1:n
        k1 = max(1,k-ml)
        k2 = min(k+mu,n)
        for j=k1:k2
            if j!=k
                fjac[k,j] = -(ONE+TWO*x[j])
            end
        end
        fjac[k,k] = TWO+FIFTN*x[k]*x[k]
    end
end

broydenbanded(n::Int) = broydentridiagonal(n)

end
