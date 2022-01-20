# calculates remainder from inversion formula for x < -1
function li_neg_rest(n::Integer, x::Float64)::Float64
    l = log(-x)
    sum = 0.0

    for r in 1:(n÷2)
        sum += l^(n - 2*r)/factorial(n - 2*r)*(2.0^(1 - 2*r) - 1.0)*zeta(2*r)
    end

    2*sum - l^n/factorial(n)
end

# calculates remainder from inversion formula for x > 1
function li_pos_rest(n::Integer, x::Float64)::Float64
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

    x == 1.0 && return zeta(n)
    x == -1.0 && return (2.0^(1 - n) - 1.0)*zeta(n)

    # transformation of x to [-1,1]
    (y, rest, sgn) = if x < -1.0
        (inv(x), li_neg_rest(n, x), isodd(n) ? 1.0 : -1.0)
    elseif x < 1.0
        (x, 0.0, 1.0)
    else # x > 1.0
        (inv(x), li_pos_rest(n, x), isodd(n) ? 1.0 : -1.0)
    end

    rest + sgn*li_series(n, y)
end
