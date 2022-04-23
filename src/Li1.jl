"""
    li1(x::Real)

Returns the real 1st order polylogarithm
``\\Re[\\operatorname{Li}_1(x)]`` of a real number ``x`` of type
`Real`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li1(0.5)
0.6931471805599453
```
"""
function li1(x::Real)
    if x < one(x)
        -log(one(x) - x)
    elseif x == one(x)
        Inf
    else # x > 1
        -log(x - one(x))
    end
end

"""
    li1(z::Complex)

Returns the complex 1st order polylogarithm
``\\operatorname{Li}_1(z)`` of a complex number ``z`` of type
`Complex`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li1(1.0 + 1.0im)
-0.0 + 1.5707963267948966im
```
"""
function li1(z::Complex)
    -clog(one(z) - z)
end
