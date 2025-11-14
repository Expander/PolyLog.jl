# returns n-th harmonic number, n > 0
function harmonic(n::Integer, ::Type{T})::T where T
    if n <= 0
        throw(DomainError(n, "harmonic not implemented for n <= 0"))
    elseif n == 1
        one(T)
    elseif issimplefloat(T) && n >= 20
        convert(T, Base.MathConstants.eulergamma + digamma(n + 1))
    else
        one(T) + sum(x -> inv(T(x)), 2:n)
    end
end
