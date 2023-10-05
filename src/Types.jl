# returns true if type of x is of a simple float type, otherwise false
issimplefloat(x::Real) = issimplefloat(typeof(x))

# returns true if T is Float16, Float32 or Float64, otherwise false
issimplefloat(::Type{T}) where T = (T in (Float16, Float32, Float64))
