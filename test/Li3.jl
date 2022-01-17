@testset "li3" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Li3.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        z        = row[1] + row[2]*1im
        expected = row[3] + row[4]*1im

        if imag(z) == 0.0
            @test PolyLog.li3(real(z)) ≈ real(expected) atol=1e-14
        end

        @test PolyLog.li3(z) ≈ expected atol=1e-14
    end

    @test PolyLog.li3(0.0) == 0.0
    @test PolyLog.li3(0.5) == 0.53721319360804020 # (-2*Pi^2*log(2) + 4*log(2)^3 + 21*zeta(3))/24
    @test PolyLog.li3(1.0) == 1.2020569031595943 # zeta(3)
    @test PolyLog.li3(-1.0) == -0.90154267736969571 # -3/4*zeta(3)
end
