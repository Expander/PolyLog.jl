# returns r.h.s. of inversion formula for x < -1:
#
# Li(n,-x) + (-1)^n Li(n,-1/x)
#    = -log(n,x)^n/n! + 2 sum(r=1:(n÷2), log(x)^(n-2r)/(n-2r)! Li(2r,-1))
function li_rem_neg(n::Integer, x::Float64)::Float64
    l = log(-x)
    l2 = l*l
    sum = 0.0

    if iseven(n)
        p = 1.0 # collects l^(2u)
        for u in 0:(n÷2 - 1)
            old_sum = sum
            sum += p*inv_fac(2*u)*neg_eta(n - 2*u)
            sum == old_sum && break
            p *= l2
        end
    else
        p = l # collects l^(2u + 1)
        for u in 0:((n - 3)÷2)
            old_sum = sum
            sum += p*inv_fac(2*u + 1)*neg_eta(n - 1 - 2*u)
            sum == old_sum && break
            p *= l2
        end
    end

    2*sum - p*inv_fac(n)
end

# returns r.h.s. of inversion formula for x > 1;
# same expression as in li_rem_neg(n,x), but with
# complex logarithm log(Complex(-x))
function li_rem_pos(n::Integer, x::Float64)::Float64
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
            sum += p*co*inv_fac(2*u)*neg_eta(n - 2*u)
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
            sum += p*co*inv_fac(2*u + 1)*neg_eta(n - 1 - 2*u)
            sum == old_sum && break
            p *= l2
            co, si = co*c2 - si*s2, si*c2 + co*s2
        end
    end

    2*sum - p*co*inv_fac(n)
end

# returns Li(n,x) using the series expansion of Li(n,x) for n > 0 and
# x ~ 1 where 0 < x < 1:
#
# Li(n,x) = sum(j=0:Inf, zeta(n-j) log(x)^j/j!)
#
# where
#
# zeta(1) = -log(-log(x)) + harmonic(n - 1)
#
# harmonic(n) = sum(k=1:n, 1/k)
function li_series_unity_pos(n::Integer, x::Float64)::Float64
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

# returns Li(n,x) using the series expansion of Li(n,x) for n < 0 and
# x ~ 1
#
# Li(n,x) = gamma(1-n) (-ln(x))^(n-1)
#           + sum(k=0:Inf, zeta(n-k) ln(x)^k/k!)
function li_series_unity_neg(n::Integer, z::ComplexF64)::ComplexF64
    clog(z) = 0.5*log(abs2(z)) + angle(z)*1.0im

    lnz = clog(z)
    lnz2 = lnz*lnz
    sum = fac(-n)*(-lnz)^(n - 1)

    if iseven(n)
        kmin = 1
        lnzk = lnz
    else
        kmin = 2
        lnzk = lnz2
        sum += zeta(n)
    end

    for k in kmin:2:typemax(n)
        term = zeta(n - k)*inv_fac(k)*lnzk
        !isfinite(term) && break
        sum_old = sum
        sum += term
        sum == sum_old && break
        lnzk *= lnz2
    end

    sum
end

# returns Li(n,x) using the naive series expansion of Li(n,x)
# for |x| < 1:
#
# Li(n,x) = sum(k=1:Inf, x^k/k^n)
function li_series_naive(n::Integer, x::Float64)::Float64
    sum = x
    xn = x*x

    for k in 2:typemax(n)
        term = xn/Float64(k)^n
        !isfinite(term) && break
        old_sum = sum
        sum += term
        sum == old_sum && break
        xn *= x
    end

    sum
end

# returns |ln(x)|^2 for all x
function ln_sqr(x::Float64)::Float64
    if x < 0.0
        log(-x)^2 + pi^2
    elseif x == 0.0
        NaN
    else
        log(x)^2
    end
end

"""
    li(n::Integer, x::Float64)::Float64

Returns the real n-th order polylogarithm
``\\Re[\\operatorname{Li}_n(x)]`` of a real number ``x`` of type
`Float64` for all integers ``n``.

The implementation for ``n < 0`` is an adaption of
[[arxiv:2010.09860](https://arxiv.org/abs/2010.09860)].

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li(10, 1.0)
1.0009945751278182
```
"""
function li(n::Integer, x::Float64)::Float64
    isnan(x) && return NaN
    isinf(x) && return -Inf
    x == 0.0 && return 0.0
    x == 1.0 && return zeta(n)
    x == -1.0 && return neg_eta(n)

    if n < 0
        # arXiv:2010.09860
        l2 = ln_sqr(x)
        if 4*pi^2*x*x < l2
            li_series_naive(n, x)
        elseif l2 < 0.512*0.512*4*pi^2
            real(li_series_unity_neg(n, Complex(x)))
        else
            (isodd(n) ? 1.0 : -1.0)*li_series_naive(n, inv(x))
        end
    elseif n == 0
        li0(x)
    elseif n == 1
        li1(x)
    elseif n == 2
        li2(x)
    elseif n == 3
        li3(x)
    elseif n == 4
        li4(x)
    else # n > 4
        # transform x to [-1,1]
        (x, rest, sgn) = if x < -1.0
            (inv(x), li_rem_neg(n, x), isodd(n) ? 1.0 : -1.0)
        elseif x < 1.0
            (x, 0.0, 1.0)
        else # x > 1.0
            (inv(x), li_rem_pos(n, x), isodd(n) ? 1.0 : -1.0)
        end

        li = if n < 20 && x > 0.75
            li_series_unity_pos(n, x)
        else
            li_series_naive(n, x)
        end

        rest + sgn*li
    end
end
