function is_equal(a::Float64, b::Float64, eps::Float64)
    return abs(a - b) <= (1.0 + max(abs(a), abs(b))) * eps
end


function is_equal(a::ComplexF64, b::ComplexF64, eps::Float64)
    return is_equal(real(a), real(b), eps) && is_equal(imag(a), imag(b), eps)
end
