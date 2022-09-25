import SpecialFunctions

"""
    relihalf(n::Integer, x::Real)

Returns the real n/2-th order polylogarithm
``\\Re[\\operatorname{Li}_{n/2}(x)]`` of a real number ``x`` of type
`Real` for all integers ``n``.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> relihalf(1, 0.5)
0.80612672304285226
```
"""
relihalf(n::Integer, x::Real) = _relihalf(n, float(x))

_relihalf(n::Integer, x::Float16) = oftype(x, _relihalf(n, Float32(x)))

_relihalf(n::Integer, x::Float32) = oftype(x, _relihalf(n, Float64(x)))

function _relihalf(n::Integer, x::Float64)::Float64
    isnan(x) && return NaN
    isinf(x) && return -Inf
    iseven(n) && return reli(n÷2, x)
    x == 0.0 && return 0.0
    x == 1.0 && n == 1 && return Inf
    x == 1.0 && return zetahalf(n)

    if n < 0
        throw(DomainError(n, "relihalf not implemented for n < 0"))
    else # n > 0
        if abs(x) < 0.75
            lihalf_series_naive(n, x)
        else
            real(lihalf(n, Complex(x)))
        end
    end
end

"""
    lihalf(n::Integer, z::Complex)

Returns the n/2-th order polylogarithm
``\\operatorname{Li}_{n/2}(x)`` of a complex number ``z`` of type
`Complex` for all integers ``n``.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> lihalf(1, 0.5)
-1.5466407024391607 - 2.7835446536610238im
```
"""
lihalf(n::Integer, z::Complex) = _lihalf(n, float(z))

lihalf(n::Integer, z::Real) = lihalf(n, Complex(z))

_lihalf(n::Integer, z::ComplexF16) = oftype(z, _lihalf(n, ComplexF32(z)))

_lihalf(n::Integer, z::ComplexF32) = oftype(z, _lihalf(n, ComplexF64(z)))

function _lihalf(n::Integer, z::ComplexF64)::ComplexF64
    isnan(z) && return NaN
    isinf(z) && return -Inf
    iseven(n) && return li(n÷2, z)
    z == 0.0 && return 0.0
    z == 1.0 && n == 1 && return Inf
    z == 1.0 && return zetahalf(n)

    if n < 0
        throw(DomainError(n, "lihalf not implemented for n < 0"))
    else # n > 0
        if abs2(z) < 0.75^2
            lihalf_series_naive(n, z)
        elseif abs2(z) >= 1.4^2
            -Complex(-1.0, 0.0)^(n/2)*lihalf_series_naive(n, inv(z)) + lihalf_rem(n, z)
        else
            lihalf_series_unity_pos(n, z)
        end
    end
end

# returns Li(n/2,z) using the naive series ezpansion of Li(n/2,z)
# for |z| < 1:
#
# Li(n/2,z) = sum(k=1:Inf, z^k/k^(n/2))
function lihalf_series_naive(n::Integer, z::ComplexOrReal)
    sum = z
    zn = z*z

    for k in 2:typemax(n)
        term = zn/Float64(k)^(n/2) # TODO: optimize
        !isfinite(term) && break
        old_sum = sum
        sum += term
        sum == old_sum && break
        zn *= z
    end

    sum
end

# convert -0.0 to 0.0
function posfp0(x::Real)
    x == zero(x) ? zero(x) : x
end

# convert -0.0 to 0.0
function posfp0(z::Complex)
    Complex(posfp0(real(z)), posfp0(imag(z)))
end

# returns Li(n/2,z) using the series expansion of Li(n/2,z) for n > 0
# and z ~ 1:
#
# Li(n/2,z) = gamma(1 - n/2) (-log(x))^(n/2-1) + sum(j=0:Inf, zeta(n/2-j) log(z)^j/j!)
function lihalf_series_unity_pos(n::Integer, z::Complex)
    l = clog(z)
    sum = gammahalf(2 - n)*(posfp0(-l))^(n/2 - 1) + zetahalf(n) # TODO: optimize
    p = one(l) # collects l^j/j!

    for j in 1:typemax(n)
        p *= l/j
        term = zetahalf(n - 2*j)*p
        !isfinite(term) && break
        old_sum = sum;
        sum += term
        sum == old_sum && break
    end

    sum
end

# returns r.h.s. of inversion formula for complex z:
#
# Li(s,z) + (-1)^s Li(s,1/z)
#    = (2pi*i)^s/gamma(s)*zeta(1 - s, 1/2 + log(-z)/(2pi*i))
function lihalf_rem(n::Integer, z::Complex)
    s = n/2
    tpi = 2*pi*1.0im
    tpi^s/gammahalf(n)*SpecialFunctions.zeta(1 - s, 1/2 + log(posfp0(-z))/tpi)
end
