# converts -0.0 to 0.0
function convert_minus_0(z)
    Complex(
        real(z) == -0.0 ? 0.0 : real(z),
        imag(z) == -0.0 ? 0.0 : imag(z)
    )
end

# returns complex logarithm of z;
# handles case when imag(z) == -0.0
function clog(z::Complex)
    if iszero(imag(z)) && real(z) > zero(real(z))
        Complex(log(real(z)), zero(real(z)))
    elseif iszero(imag(z))
        Complex(log(real(-z)), convert(typeof(real(z)),pi))
    else
        log(z)
    end
end

# returns logarithm of x
function clog(x::Real)
    log(x)
end

# returns complex logarithm of 1+z;
# handles case when imag(z) == -0.0
function clog1p(z::Complex)
    if iszero(imag(z)) && real(z) > -one(real(z))
        Complex(log1p(real(z)), zero(real(z)))
    else
        log1p(convert_minus_0(z))
    end
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
