# returns n-th harmonic number, n > 0
function harmonic(n::Integer, T=Float64)
    if n <= 0
        throw(DomainError(n, "harmonic not implemented for n <= 0"))
    elseif n < 20
        sum = one(T)
        for k in 2:n
            sum += inv(k)
        end
        sum
    else
        eulergamma + digamma(n + 1)
    end
end
