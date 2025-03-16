# returns n-th harmonic number, n > 0
function harmonic(n::Integer, ::Type{T})::T where T
    if issimplefloat(T)
        convert(T, harmonic_f64(n))
    else
        harmonic_big(n)
    end
end

function harmonic_f64(n::Integer)::Float64
    if n <= 0
        throw(DomainError(n, "harmonic not implemented for n <= 0"))
    elseif n == 1
        one(Float64)
    elseif n < 20
        one(Float64) + sum(inv, 2:n)
    else
        EULERGAMMA_F64 + digamma(n + 1)
    end
end

function harmonic_big(n::Integer)::BigFloat
    if n <= 0
        throw(DomainError(n, "harmonic not implemented for n <= 0"))
    elseif n == 1
        one(BigFloat)
    else
        one(BigFloat) + sum(x -> inv(big(x)), 2:n)
    end
end
