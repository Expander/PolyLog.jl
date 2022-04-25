# return entries from complex n x 2 matrix for which first element is real
function filter_real(data)
    z1 = [real(data[i,1]) for i in 1:size(data, 1) if imag(data[i,1]) == 0.0]
    z2 = [real(data[i,2]) for i in 1:size(data, 1) if imag(data[i,1]) == 0.0]
    hcat(z1, z2)
end

# return entries from complex n x 2 matrix for which |Im(z)| larger than smallest value
function filter_ComplexF16(data)
    z1 = [data[i,1] for i in 1:size(data, 1) if abs(imag(data[i,1])) > 0.5e-7]
    z2 = [data[i,2] for i in 1:size(data, 1) if abs(imag(data[i,1])) > 0.5e-7]
    hcat(z1, z2)
end

# apply function row-wise to first element and compare with second element
function test_function_on_data(fn, data, atol, rtol)
    for i in 1:size(data, 1)
        @test fn(data[i,1]) ≈ data[i,2] atol=atol rtol=rtol
    end
end
