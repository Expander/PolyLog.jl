@testset "li0" begin
    data = read_from(joinpath(@__DIR__, "data", "Li0.txt"))

    for d in eachrow(data)
        z = d[1]
        expected = d[2]

        if imag(z) == 0.0
            @test PolyLog.li0(real(z)) ≈ real(expected) atol=1e-14
        end

        @test PolyLog.li0(z) ≈ expected atol=1e-14
    end

    @test PolyLog.li0(1.0) == Inf
    @test isnan(PolyLog.li0(1.0 + 0.0im))
end
