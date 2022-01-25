# Table[PolyLog[n,-1], {n,1,54}]
const li_minus_1_coeff = (
    -0.69314718055994531, -0.82246703342411322, -0.90154267736969571,
    -0.94703282949724592, -0.97211977044690931, -0.98555109129743510,
    -0.99259381992283028, -0.99623300185264790, -0.99809429754160533,
    -0.99903950759827157, -0.99951714349806075, -0.99975768514385819,
    -0.99987854276326512, -0.99993917034597972, -0.99996955121309924,
    -0.99998476421490611, -0.99999237829204101, -0.99999618786961011,
    -0.99999809350817168, -0.99999904661158152, -0.99999952325821554,
    -0.99999976161323082, -0.99999988080131844, -0.99999994039889239,
    -0.99999997019885696, -0.99999998509923200, -0.99999999254955048,
    -0.99999999627475340, -0.99999999813736942, -0.99999999906868228,
    -0.9999999995343403 , -0.9999999997671699 , -0.9999999998835849 ,
    -0.9999999999417924 , -0.9999999999708962 , -0.9999999999854481 ,
    -0.9999999999927240 , -0.9999999999963620 , -0.9999999999981810 ,
    -0.9999999999990905 , -0.9999999999995453 , -0.9999999999997726 ,
    -0.9999999999998863 , -0.9999999999999432 , -0.9999999999999716 ,
    -0.9999999999999858 , -0.9999999999999929 , -0.9999999999999964 ,
    -0.9999999999999982 , -0.9999999999999991 , -0.9999999999999996 ,
    -0.9999999999999998 , -0.9999999999999999 , -0.9999999999999999
)

# returns Li(n,-1) = (2.0^(1 - n) - 1.0)*zeta(n) for n > 0
function li_minus_1(n::Integer)::Float64
    if n < 1
        throw(DomainError(n, "li_minus_1 not implemented for n < 1"))
    elseif n <= length(li_minus_1_coeff)
        li_minus_1_coeff[n]
    else
        -1.0
    end
end

# calculates remainder from inversion formula for x < -1
function li_neg_rest(n::Integer, x::Float64)::Float64
    l = log(-x)
    l2 = l*l
    sum = 0.0

    if iseven(n)
        p = 1.0 # collects l^(2u)
        for u in 0:(n÷2 - 1)
            old_sum = sum
            sum += p*inverse_factorial(2*u)*li_minus_1(n - 2*u)
            sum == old_sum && break
            p *= l2
        end
    else
        p = l # collects l^(2u + 1)
        for u in 0:((n - 3)÷2)
            old_sum = sum
            sum += p*inverse_factorial(2*u + 1)*li_minus_1(n - 1 - 2*u)
            sum == old_sum && break
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
        s2, c2 = sincos(2.0*arg)
        co = 1.0 # collects cos(2*u*arg)
        si = 0.0 # collects sin(2*u*arg)
        for u in 0:(n÷2 - 1)
            old_sum = sum
            sum += p*co*inverse_factorial(2*u)*li_minus_1(n - 2*u)
            sum == old_sum && break
            p *= l2
            co, si = co*c2 - si*s2, si*c2 + co*s2
        end
    else
        p = mag # collects mag^(2u + 1)
        s, c = sincos(arg)
        s2, c2 = 2.0*s*c, 2.0*c*c - 1.0 # sincos(2*arg)
        co = c # collects cos((2*u + 1)*arg)
        si = s # collects sin((2*u + 1)*arg)
        for u in 0:((n - 3)÷2)
            old_sum = sum
            sum += p*co*inverse_factorial(2*u + 1)*li_minus_1(n - 1 - 2*u)
            sum == old_sum && break
            p *= l2
            co, si = co*c2 - si*s2, si*c2 + co*s2
        end
    end

    2*sum - p*co*inverse_factorial(n)
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
`Float64` for integer ``n\\geq 0``.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li(10, 1.0)
1.0009945751278182
```
"""
function li(n::Integer, x::Float64)::Float64
    n < 0 && throw(DomainError(n, "li(n,x) not implemented for n < 0"))
    n == 0 && return li0(x)
    n == 1 && return li1(x)
    n == 2 && return li2(x)
    n == 3 && return li3(x)
    n == 4 && return li4(x)
    x == 1.0 && return zeta(n)
    x == -1.0 && return li_minus_1(n)
    isnan(x) && return NaN

    # transform x to [-1,1]
    (x, rest, sgn) = if x < -1.0
        (inv(x), li_neg_rest(n, x), isodd(n) ? 1.0 : -1.0)
    elseif x < 1.0
        (x, 0.0, 1.0)
    else # x > 1.0
        (inv(x), li_pos_rest(n, x), isodd(n) ? 1.0 : -1.0)
    end

    li = if n < 20 && x > 0.75
        li_series_one(n, x)
    else
        li_series(n, x)
    end

    rest + sgn*li
end
