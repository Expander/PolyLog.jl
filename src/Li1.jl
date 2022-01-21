"""
    li1(x::Float64)::Float64

Returns the real 1st order polylogarithm
``\\Re[\\operatorname{Li}_1(x)]`` of a real number ``x`` of type
`Float64`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li1(0.5)
0.6931471805599453
```
"""
function li1(x::Float64)::Float64
    if x < 1.0
        -log(1.0 - x)
    elseif x == 1.0
        Inf
    else # x > 1.0
        -log(x - 1.0)
    end
end
