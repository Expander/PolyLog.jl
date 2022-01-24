@testset "li1" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Li1.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        z        = row[1] + row[2]*1im
        expected = row[3] + row[4]*1im

        if imag(z) == 0.0
            @test PolyLog.li1(real(z)) ≈ real(expected) atol=1e-14
        end

        @test PolyLog.li1(z) ≈ expected atol=1e-14
    end

    @test PolyLog.li1(1.0) == Inf
end
