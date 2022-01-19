function li(n::Integer, x::Float64)::Float64
    if n < 1
        throw(DomainError(n, "li(n,x) undefined for n < 1"))
    end

    0.0
end
