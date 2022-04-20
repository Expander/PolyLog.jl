@testset "li4" begin
    data = read_from(joinpath(@__DIR__, "data", "Li4.txt"))

    for d in eachrow(data)
        z = d[1]
        expected = d[2]

        if imag(z) == 0.0
            @test PolyLog.li4(real(z)) ≈ real(expected) atol=1e-14
        end

        @test PolyLog.li4(z) ≈ expected atol=1e-14
    end
end
