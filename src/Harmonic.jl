# returns n-th harmonic number, n > 0
function harmonic(n::Integer)::Float64
    if n <= 0
        throw(DomainError(n, "harmonic not implemented for n <= 0"))
    elseif n < 20
        sum = 1.0
        for k in 2:n
            sum += 1/k
        end
        sum
    else
        eulergamma + digamma(n + 1)
    end
end
