"""
    li0(x::Float64)::Float64

Returns the real 0th order polylogarithm
``\\Re[\\operatorname{Li}_0(x)]`` of a real number ``x`` of type
`Float64`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li0(0.5)
1.0
```
"""
function li0(x::Float64)::Float64
    x/(1 - x)
end
