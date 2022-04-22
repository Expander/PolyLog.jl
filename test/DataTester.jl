function test_cmpl_function_on_data(fn, filename, atol)
    data = read_from(filename)
    for i in 1:size(data, 1)
        @test fn(data[i,1]) ≈ data[i,2] atol=atol
    end
end

function test_real_function_on_data(fn, filename, atol)
    data = read_from(filename)
    for i in 1:size(data, 1)
        z = data[i,1]
        if imag(z) == 0.0
            @test fn(real(z)) ≈ real(data[i,2]) atol=atol
        end
    end
end
