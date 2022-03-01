const ComplexOrReal{T} = Union{T,Complex{T}}

# returns r.h.s. of inversion formula for complex z:
#
# Li(n,-z) + (-1)^n Li(n,-1/z)
#    = -log(n,z)^n/n! + 2 sum(r=1:(n÷2), log(z)^(n-2r)/(n-2r)! Li(2r,-1))
function li_rem(n::Integer, z::ComplexF64)::ComplexF64
    l = log(-z)
    l2 = l*l;
    kmax = iseven(n) ? n÷2 : (n - 1)÷2
    p = iseven(n) ? 1.0 + 0.0im : l
    sum = 0.0 + 0.0im

    for k in kmax:-1:1
        ifac = inv_fac(n - 2*k)
        ifac == 0.0 && return 2*sum
        old_sum = sum
        sum += neg_eta(2*k)*ifac*p
        sum == old_sum && break
        p *= l2
    end

    2*sum - p*inv_fac(n)
end

# returns r.h.s. of inversion formula for real x < -1:
#
# Li(n,-x) + (-1)^n Li(n,-1/x)
#    = -log(n,x)^n/n! + 2 sum(r=1:(n÷2), log(x)^(n-2r)/(n-2r)! Li(2r,-1))
function li_rem(n::Integer, x::Float64)::Float64
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

# returns Li(n,z) using the series expansion of Li(n,z) for n > 0 and
# z ~ 1:
#
# Li(n,z) = sum(j=0:Inf, zeta(n-j) log(z)^j/j!)
#
# where
#
# zeta(1) = -log(-log(z)) + harmonic(n - 1)
#
# harmonic(n) = sum(k=1:n, 1/k)
function li_series_unity_pos(n::Integer, z::ComplexOrReal)
    l = log(z)
    sum = zeta(n)
    p = 1.0 # collects l^j/j!

    for j in 1:(n - 2)
        p *= l/j
        old_sum = sum;
        sum += zeta(n - j)*p
        sum == old_sum && break
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

# returns Li(n,z) using the series expansion of Li(n,z) for n < 0 and
# z ~ 1
#
# Li(n,z) = gamma(1-n) (-ln(z))^(n-1)
#           + sum(k=0:Inf, zeta(n-k) ln(z)^k/k!)
function li_series_unity_neg(n::Integer, z::ComplexF64)::ComplexF64
    lnz = log(z)
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
function li_series_naive(n::Integer, z::ComplexOrReal)
    sum = z
    zn = z*z

    for k in 2:typemax(n)
        term = zn/Float64(k)^n
        !isfinite(term) && break
        old_sum = sum
        sum += term
        sum == old_sum && break
        zn *= z
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

# returns -(-1)^n
oddsgn(n) = isodd(n) ? 1.0 : -1.0

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
            oddsgn(n)*li_series_naive(n, inv(x))
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
            (inv(x), li_rem(n, x), oddsgn(n))
        elseif x < 1.0
            (x, 0.0, 1.0)
        else # x > 1.0
            (inv(x), real(li_rem(n, Complex(x))), oddsgn(n))
        end

        li = if n < 20 && x > 0.75
            li_series_unity_pos(n, x)
        else
            li_series_naive(n, x)
        end

        rest + sgn*li
    end
end

"""
    li(n::Integer, x::ComplexF64)::ComplexF64

Returns the complex n-th order polylogarithm
``\\operatorname{Li}_n(x)`` of a complex number ``z`` of type
`ComplexF64` for all integers ``n``.

The implementation for ``n < 0`` is an adaption of
[[arxiv:2010.09860](https://arxiv.org/abs/2010.09860)].

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li(10, Complex(1.0, 1.0))
1.0009945751278182
```
"""
function li(n::Integer, z::ComplexF64)::ComplexF64
    isnan(z) && return NaN
    isinf(z) && return -Inf
    z == 0.0 && return 0.0
    z == 1.0 && return zeta(n)
    z == -1.0 && return neg_eta(n)

    if n < 0
        l2 = abs2(log(z))
        if 4*pi^2*abs2(z) < l2
            li_series_naive(n, z)
        elseif l2 < 0.512*0.512*4*pi^2
            li_series_unity_neg(n, z)
        else
            sqrtz = sqrt(z)
            2.0^(n - 1)*(li(n, sqrtz) + li(n, -sqrtz))
        end
    elseif n == 0
        li0(z)
    elseif n == 1
        li1(z)
    elseif n == 2
        li2(z)
    elseif n == 3
        li3(z)
    elseif n == 4
        li4(z)
    elseif n == 5
        li5(z)
    elseif n == 6
        li6(z)
    elseif abs2(z) <= 0.75*0.75
        li_series_naive(n, z)
    elseif abs2(z) >= 1.4*1.4
        oddsgn(n)*li_series_naive(n, 1.0/z) + li_rem(n, z)
    else
        li_series_unity_pos(n, z)
    end
end
