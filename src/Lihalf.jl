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
            return lihalf_series_naive(n, x)
        elseif ln_sqr(x) < (2*pi)^2
            return real(lihalf_series_unity_pos(n, Complex(x)))
        end
        throw(DomainError(n, "relihalf not implemented for n > 0 and abs(x) >= 0.75"))
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

# returns Li(n/2,z) using the series expansion of Li(n/2,z) for n > 0
# and z ~ 1:
#
# Li(n/2,z) = gamma(1 - n/2) (-log(x))^(n/2-1) + sum(j=0:Inf, zeta(n/2-j) log(z)^j/j!)
function lihalf_series_unity_pos(n::Integer, z::Complex)
    l = clog(z)
    sum = gammahalf(2 - n)*(-l)^((n - 1)÷2)/sqrt(-l) + zetahalf(n)
    p = 1.0 # collects l^j/j!

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
