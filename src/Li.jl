# calculates remainder from inversion formula
function li_rest(n::Integer, x::Float64)::Float64
    0.0
end

# naive series expansion of Li(n,x) for |x| < 1
function li_series(n::Integer, x::Float64)::Float64
    sum = x

    for k in 2:typemax(n)
        old_sum = sum
        sum += x^k/Float64(k)^n
        sum == old_sum && break
    end

    sum
end

function li(n::Integer, x::Float64)::Float64
    if n < 1
        throw(DomainError(n, "li(n,x) undefined for n < 1"))
    end

    # transformation on [-1,1]
    (y, rest, sgn) = if (abs(x) <= 1.0)
        (x, 0.0, 1.0)
    else # abs(x) > 1
        (inv(x), li_rest(n, x), iseven(n - 1) ? 1.0 : -1.0)
    end

    li_series(n, y)
end
