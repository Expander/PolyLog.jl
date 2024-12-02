"""
    li0(z::Number)

Returns the 0th order polylogarithm ``\\operatorname{Li}_0(z) = z/(1-z)``
of a number ``z`` of type `Number`.


Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li0(0.25)
0.3333333333333333

julia> li0(BigFloat("0.25"))
0.3333333333333333333333333333333333333333333333333333333333333333333333333333348
```
"""
function li0(z::Number)
    if iszero(z)
        z
    else
        z/(1 - z)
    end
end
