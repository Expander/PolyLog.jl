@testset "li2" begin
    data = read_from(joinpath(@__DIR__, "data", "Li2.txt"))

    for i in 1:size(data, 1)
        z = data[i,1]
        expected = data[i,2]

        if imag(z) == 0.0
            @test PolyLog.li2(real(z)) ≈ real(expected) atol=1e-14
        end

        @test PolyLog.li2(z) ≈ expected atol=1e-14
    end

    @test PolyLog.li2(1.0) == 1.6449340668482264
    @test real(PolyLog.li2(-1.08371e-08 + 1.32716e-24*im)) == -1.0837099970639316e-08
    @test imag(PolyLog.li2(-1.08371e-08 + 1.32716e-24*im)) ==  1.3271599928087172e-24
end
