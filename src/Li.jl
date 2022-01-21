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
    l = log(Complex(-x))
    sum = 0.0

    for r in 1:(n÷2)
        sum += l^(n - 2*r)/factorial(n - 2*r)*(2.0^(1 - 2*r) - 1.0)*zeta(2*r)
    end

    real(2*sum - l^n/factorial(n))
end

function harmonic(n::Integer)
    sum = 1.0

    for k in 2:n
        sum += 1/k
    end

    sum
end

# series expansion of Li(n,x) for x ~ 1, x < 1
function li_series_one(n::Integer, x::Float64)::Float64
    l = log(x)
    sum = zeta(n)
    p = 1.0 # collects l^j/j!

    for j in 1:(n + 1)
        p *= l/j
        if j == n - 1
            sum += (harmonic(n - 1) - log(-l))*p
        else
            sum += zeta(n - j)*p
        end
    end

    l2 = l*l

    for j in (n + 3):2:typemax(n)
        p *= l2/((j - 1)*j)
        old_sum = sum
        sum += zeta(n - j)*p
        sum == old_sum && break
    end

    sum
end

# naive series expansion of Li(n,x) for |x| < 1
function li_series(n::Integer, x::Float64)::Float64
    sum = x
    xn = x*x

    for k in 2:typemax(n)
        old_sum = sum
        sum += xn/Float64(k)^n
        sum == old_sum && break
        xn *= x
    end

    sum
end

function li(n::Integer, x::Float64)::Float64
    if n < 0
        throw(DomainError(n, "li(n,x) not implemented for n < 0"))
    end

    n == 0 && return x/(1 - x)
    n == 1 && return -log(1.0 - x)
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

    li = x > 0.75 ? li_series_one(n, y) : li_series(n, y)

    rest + sgn*li
end
