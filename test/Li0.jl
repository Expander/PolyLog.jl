@testset "li0" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Li0.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        z        = row[1] + row[2]*1im
        expected = row[3] + row[4]*1im

        if imag(z) == 0.0
            @test PolyLog.li0(real(z)) ≈ real(expected) atol=1e-14
        end

        @test PolyLog.li0(z) ≈ expected atol=1e-14
    end

    @test PolyLog.li0(1.0) == Inf
    @test isnan(PolyLog.li0(1.0 + 0.0im))
end
