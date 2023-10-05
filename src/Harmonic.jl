# returns n-th harmonic number, n > 0
function harmonic(n::Integer, T)
    if issimplefloat(T)
        convert(T, harmonic_f64(n))
    else
        harmonic_big(n)
    end
end

function harmonic_f64(n::Integer)::Float64
    if n <= 0
        throw(DomainError(n, "harmonic not implemented for n <= 0"))
    elseif n < 20
        sum = one(Float64)
        for k in 2:n
            sum += inv(k)
        end
        sum
    else
        eulergamma + digamma(n + 1)
    end
end

function harmonic_big(n::Integer)::BigFloat
    if n <= 0
        throw(DomainError(n, "harmonic not implemented for n <= 0"))
    else
        sum = one(BigFloat)
        for k in 2:n
            sum += inv(big(k))
        end
        sum
    end
end
