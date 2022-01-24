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

"""
    li0(x::ComplexF64)::ComplexF64

Returns the complex 0th order polylogarithm
``\\operatorname{Li}_0(z)`` of a complex number ``z`` of type
`ComplexF64`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li0(1.0 + 1.0im)
-1.0 + 1.0im
```
"""
function li0(z::ComplexF64)::ComplexF64
    z/(1 - z)
end
