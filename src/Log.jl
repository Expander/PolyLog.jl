# returns complex logarithm of z;
# handles case when imag(z) == -0.0
function clog(z::Complex)
    ang = angle(z)
    Complex(log(abs(z)), iszero(imag(z)) && ang < zero(z.re) ? -ang : ang)
end

# returns logarithm of x
function clog(x::Real)
    log(x)
end

# returns complex logarithm of 1+z;
# handles case when imag(z) == -0.0
function clog1p(z::Complex)
    clog(one(z) + z)
end

# returns |ln(x)|^2 for all x
function ln_sqr(x::Real)::Real
    if x < zero(x)
        abs2(log(complex(x)))
    elseif iszero(x)
        oftype(x, Inf)
    else
        log(x)^2
    end
end

# returns |ln(x)|^2 for all x
function ln_sqr(x::Float64)::Float64
    if x < zero(x)
        log(-x)^2 + pi^2
    elseif iszero(x)
        Inf
    else
        log(x)^2
    end
end
