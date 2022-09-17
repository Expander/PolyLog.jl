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
    x == 0.0 && return 0.0

    if n < 0
        throw(DomainError(n, "relihalf not implemented for n < 0"))
    elseif n == 0
        li0(x)
    else
        throw(DomainError(n, "relihalf not implemented for n > 0"))
    end
end
