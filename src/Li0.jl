"""
    li0(z::Number)

Returns the 0th order polylogarithm ``\\operatorname{Li}_0(z) = z/(1-z)``
of a number ``z`` of type `Number`.


Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li0(0.5)
1.0
```
"""
li0(z::Number) = z/(1 - z)
