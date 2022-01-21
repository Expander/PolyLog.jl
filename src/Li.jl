# returns Li(n,-1)
li_minus_1(n::Integer)::Float64 = (2.0^(1 - n) - 1.0)*zeta(n)

# calculates remainder from inversion formula for x < -1
function li_neg_rest(n::Integer, x::Float64)::Float64
    l = log(-x)
    l2 = l*l
    sum = 0.0

    if iseven(n)
        p = 1.0 # collects l^(2u)
        for u in 0:(n÷2 - 1)
            sum += p*inverse_factorial(2*u)*li_minus_1(n - 2*u)
            p *= l2
        end
    else
        p = l # collects l^(2u + 1)
        for u in 0:((n - 3)÷2)
            sum += p*inverse_factorial(2*u + 1)*li_minus_1(n - 1 - 2*u)
            p *= l2
        end
    end

    2*sum - p*inverse_factorial(n)
end

# calculates remainder from inversion formula for x > 1
function li_pos_rest(n::Integer, x::Float64)::Float64
    l = log(x)
    mag = hypot(l, pi) # |log(-x)|
    arg = atan(pi, l)  # angle(log(-x))
    l2 = mag*mag       # |log(-x)|^2
    sum = 0.0

    if iseven(n)
        p = 1.0 # collects mag^(2u)
        for u in 0:(n÷2 - 1)
            sum += p*cos(2*u*arg)*inverse_factorial(2*u)*li_minus_1(n - 2*u)
            p *= l2
        end
    else
        p = mag # collects mag^(2u + 1)
        for u in 0:((n - 3)÷2)
            sum += p*cos((2*u + 1)*arg)*inverse_factorial(2*u + 1)*li_minus_1(n - 1 - 2*u)
            p *= l2
        end
    end

    2*sum - p*cos(n*arg)*inverse_factorial(n)
end

# returns n-th harmonic number
function harmonic(n::Integer)::Float64
    sum = 1.0

    for k in 2:n
        sum += 1/k
    end

    sum
end

# series expansion of Li(n,x) for x ~ 1, 0 < x < 1
function li_series_one(n::Integer, x::Float64)::Float64
    l = log(x)
    sum = zeta(n)
    p = 1.0 # collects l^j/j!

    for j in 1:(n - 2)
        p *= l/j
        sum += zeta(n - j)*p
    end

    p *= l/(n - 1)
    sum += (harmonic(n - 1) - log(-l))*p

    p *= l/n
    sum += zeta(0)*p

    p *= l/(n + 1)
    sum += zeta(-1)*p

    l2 = l*l

    for j in (n + 3):2:typemax(n)
        p *= l2/((j - 1)*j)
        old_sum = sum
        sum += zeta(n - j)*p
        sum == old_sum && break
    end

    sum
end

# naive series expansion of Li(n,x) for |x| < 1
function li_series(n::Integer, x::Float64)::Float64
    sum = x
    xn = x*x

    for k in 2:typemax(n)
        old_sum = sum
        sum += xn/Float64(k)^n
        sum == old_sum && break
        xn *= x
    end

    sum
end

"""
    li(n::Integer, x::Float64)::Float64

Returns the real n-th order polylogarithm
``\\Re[\\operatorname{Li}_n(x)]`` of a real number ``x`` of type
`Float64`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li(10, 1.0)
1.0009945751278182
```
"""
function li(n::Integer, x::Float64)::Float64
    if n < 0
        throw(DomainError(n, "li(n,x) not implemented for n < 0"))
    end

    n == 0 && return li0(x)
    n == 1 && return li1(x)
    x == 1.0 && return zeta(n)
    x == -1.0 && return li_minus_1(n)

    # transformation of x to [-1,1]
    (y, rest, sgn) = if x < -1.0
        (inv(x), li_neg_rest(n, x), isodd(n) ? 1.0 : -1.0)
    elseif x < 1.0
        (x, 0.0, 1.0)
    else # x > 1.0
        (inv(x), li_pos_rest(n, x), isodd(n) ? 1.0 : -1.0)
    end

    li = x > 0.75 ? li_series_one(n, y) : li_series(n, y)

    rest + sgn*li
end
