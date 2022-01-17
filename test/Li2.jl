@testset "li2" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Li2.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        z        = row[1] + row[2]*1im
        expected = row[3] + row[4]*1im

        if imag(z) == 0.0
            @test PolyLog.li3(real(z)) ≈ real(expected) atol=1e-14
        end

        @test PolyLog.li2(z) ≈ expected atol=1e-14
    end

    @test PolyLog.li2(1.0) == 1.6449340668482264
    @test real(PolyLog.li2(-1.08371e-08 + 1.32716e-24*im)) == -1.0837099970639316e-08
    @test imag(PolyLog.li2(-1.08371e-08 + 1.32716e-24*im)) ==  1.3271599928087172e-24
end
