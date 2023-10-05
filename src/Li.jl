const ComplexOrReal{T} = Union{T,Complex{T}}

"""
    reli(n::Integer, x::Real)

Returns the real n-th order polylogarithm
``\\Re[\\operatorname{Li}_n(x)]`` of a real number ``x`` of type
`Real` for all integers ``n``.

The implementation for ``n < 0`` is an adaptation of
[[arxiv:2010.09860](https://arxiv.org/abs/2010.09860)].

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> reli(10, 1.0)
1.0009945751278182

julia> reli(10, BigFloat("1.0"))
1.000994575127818085337145958900319017006019531564477517257788994636291465151908
```
"""
reli(n::Integer, x::Real) = _reli(n, float(x))

_reli(n::Integer, x::Float16) = oftype(x, _reli(n, Float32(x)))

_reli(n::Integer, x::Float32) = oftype(x, _reli(n, Float64(x)))

function _reli(n::Integer, x::Real)
    isnan(x) && return oftype(x, NaN)
    isinf(x) && return oftype(x, -Inf)
    x == zero(x) && return zero(x)
    x == one(x) && return zeta(n, typeof(x))
    x == -one(x) && return neg_eta(n, typeof(x))

    if n < 0
        # arXiv:2010.09860
        l2 = ln_sqr(x)
        if (2*pi)^2*x*x < l2
            li_series_naive(n, x)
        elseif l2 < (512/1000*2*pi)^2
            real(li_series_unity_neg(n, Complex(x)))
        else
            oddsgn(n, typeof(x))*li_series_naive(n, inv(x))
        end
    elseif n == 0
        li0(x)
    elseif n == 1
        reli1(x)
    elseif n == 2
        reli2(x)
    elseif issimplefloat(x) && n == 3
        reli3(x)
    elseif issimplefloat(x) && n == 4
        reli4(x)
    else # n > 4
        # transform x to [-1,1]
        (x, rest, sgn) = if x < -one(x)
            (inv(x), reli_rem(n, x), oddsgn(n, typeof(x)))
        elseif x < one(x)
            (x, zero(x), one(x))
        else # x > one(x)
            (inv(x), real(li_rem(n, Complex(x))), oddsgn(n, typeof(x)))
        end

        li = if n < 20 && x > 3/4
            li_series_unity_pos(n, x)
        else
            li_series_naive(n, x)
        end

        rest + sgn*li
    end
end

"""
    li(n::Integer, z::Complex)

Returns the complex n-th order polylogarithm
``\\operatorname{Li}_n(z)`` of a complex number ``z`` of type
`Complex` for all integers ``n``.

The implementation for ``n < 0`` is an adaptation of
[[arxiv:2010.09860](https://arxiv.org/abs/2010.09860)].

If only real arguments ``z\\in\\mathbb{R}`` are considered and one is
interested only in the real part of the polylogarithm,
``\\Re[\\operatorname{Li}_n(z)]``, refer to the function
[`reli`](@ref), which may be a faster alternative.

See also [`reli`](@ref).

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li(10, 1.0 + 1.0im)
0.9999619510320734 + 1.0019864330842578im
```
"""
li(n::Integer, z::Complex) = _li(n, float(z))

li(n::Integer, z::Real) = li(n, Complex(z))

_li(n::Integer, z::ComplexF16) = oftype(z, _li(n, ComplexF32(z)))

_li(n::Integer, z::ComplexF32) = oftype(z, _li(n, ComplexF64(z)))

function _li(n::Integer, z::ComplexF64)::ComplexF64
    isnan(z) && return NaN
    isinf(z) && return -Inf

    if iszero(imag(z))
        if real(z) <= one(typeof(real(z))) || n <= 0
            return Complex(reli(n, real(z)))
        else
            return Complex(reli(n, real(z)), -pi*inv_fac(n - 1, Float64)*log(real(z))^(n - 1))
        end
    end

    if n < 0
        l2 = abs2(log(z))
        if 4*pi^2*abs2(z) < l2
            li_series_naive(n, z)
        elseif l2 < (512/1000*2*pi)^2
            li_series_unity_neg(n, z)
        else
            sqrtz = sqrt(z)
            exp2(n - 1)*(li(n, sqrtz) + li(n, -sqrtz))
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
    elseif abs2(z) <= (3/4)^2
        li_series_naive(n, z)
    elseif abs2(z) >= (7/5)^2
        oddsgn(n, typeof(real(z)))*li_series_naive(n, inv(z)) + li_rem(n, z)
    else
        li_series_unity_pos(n, z)
    end
end

# returns -(-1)^n of type T
function oddsgn(n, ::Type{T})::T where T
    isodd(n) ? one(T) : -one(T)
end

# returns r.h.s. of inversion formula for complex z:
#
# Li(n,-z) + (-1)^n Li(n,-1/z)
#    = -log(n,z)^n/n! + 2 sum(r=1:(n÷2), log(z)^(n-2r)/(n-2r)! Li(2r,-1))
function li_rem(n::Integer, z::Complex{T})::Complex{T} where T
    l = clog(-z)
    l2 = l*l;
    kmax = iseven(n) ? n÷2 : (n - 1)÷2
    p = iseven(n) ? one(z) : l
    sum = zero(z)

    for k in kmax:-1:1
        ifac = inv_fac(n - 2*k, T)
        iszero(ifac) && return 2*sum
        old_sum = sum
        sum += neg_eta(2*k, T)*ifac*p
        p *= l2
        sum == old_sum && break
    end

    2*sum - p*inv_fac(n, T)
end

# returns r.h.s. of inversion formula for real x < -1:
#
# Li(n,-x) + (-1)^n Li(n,-1/x)
#    = -log(n,x)^n/n! + 2 sum(r=1:(n÷2), log(x)^(n-2r)/(n-2r)! Li(2r,-1))
function reli_rem(n::Integer, x::Real)
    l = log(-x)
    l2 = l*l
    sum = zero(x)

    if iseven(n)
        p = one(x) # collects l^(2u)
        for u in 0:(n÷2 - 1)
            old_sum = sum
            sum += p*inv_fac(2*u, typeof(x))*neg_eta(n - 2*u, typeof(x))
            p *= l2
            sum == old_sum && break
        end
    else
        p = l # collects l^(2u + 1)
        for u in 0:((n - 3)÷2)
            old_sum = sum
            sum += p*inv_fac(2*u + 1, typeof(x))*neg_eta(n - 1 - 2*u, typeof(x))
            p *= l2
            sum == old_sum && break
        end
    end

    2*sum - p*inv_fac(n, typeof(x))
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
function li_series_unity_pos(n::Integer, z::ComplexOrReal{T}) where T
    l = clog(z)
    sum = zeta(n, T)
    p = one(z) # collects l^j/j!

    for j in 1:(n - 2)
        p *= l/j
        old_sum = sum;
        sum += zeta(n - j, T)*p
        sum == old_sum && break
    end

    p *= l/(n - 1)
    sum += (harmonic(n - 1, T) - clog(-l))*p

    p *= l/n
    sum += zeta(0, T)*p

    p *= l/(n + 1)
    sum += zeta(-1, T)*p

    l2 = l*l

    for j in (n + 3):2:typemax(n)
        p *= l2/((j - 1)*j)
        old_sum = sum
        sum += zeta(n - j, T)*p
        sum == old_sum && break
    end

    sum
end

# returns Li(n,z) using the series expansion of Li(n,z) for n < 0 and
# z ~ 1
#
# Li(n,z) = gamma(1-n) (-ln(z))^(n-1)
#           + sum(k=0:Inf, zeta(n-k) ln(z)^k/k!)
function li_series_unity_neg(n::Integer, z::Complex{T})::Complex{T} where T
    l = clog(z)
    l2 = l*l
    sum = fac(-n, T)*(-l)^(n - 1)

    if iseven(n)
        kmin = 1
        lk = l
    else
        kmin = 2
        lk = l2
        sum += zeta(n, T)
    end

    for k in kmin:2:typemax(n)
        term = zeta(n - k, T)*inv_fac(k, T)*lk
        !isfinite(term) && break
        sum_old = sum
        sum += term
        sum == sum_old && break
        lk *= l2
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
        term = zn/oftype(real(z), k)^n
        !isfinite(term) && break
        old_sum = sum
        sum += term
        sum == old_sum && break
        zn *= z
    end

    sum
end
