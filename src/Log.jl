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
