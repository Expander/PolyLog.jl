# returns complex logarithm of z;
# handles case when imag(z) == -0.0
function clog(z::ComplexF64)::ComplexF64
    ang::Float64 = angle(z)
    Complex(0.5*log(abs2(z)), imag(z) == 0.0 && ang < 0.0 ? -ang : ang)
end

# returns logarithm of x
function clog(x::Float64)::Float64
    log(x)
end

# returns |ln(x)|^2 for all x
function ln_sqr(x::Float64)::Float64
    if x < 0.0
        log(-x)^2 + pi^2
    elseif x == 0.0
        Inf
    else
        log(x)^2
    end
end
