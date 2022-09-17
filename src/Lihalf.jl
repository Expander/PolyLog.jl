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

    if n < 0
        throw(DomainError(n, "relihalf not implemented for n < 0"))
    else # n > 0
        if abs(x) < 0.75
            return lihalf_series_naive(n, x)
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
