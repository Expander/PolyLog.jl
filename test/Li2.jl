@testset "li2" begin
    test_cmpl_function_on_data(z -> PolyLog.li2(z), joinpath(@__DIR__, "data", "Li2.txt"), 1e-14, 1e-14)
    test_real_function_on_data(z -> PolyLog.li2(z), joinpath(@__DIR__, "data", "Li2.txt"), 1e-14, 1e-14)

    @test PolyLog.li2(1.0) == 1.6449340668482264
    @test real(PolyLog.li2(-1.08371e-08 + 1.32716e-24*im)) == -1.0837099970639316e-08
    @test imag(PolyLog.li2(-1.08371e-08 + 1.32716e-24*im)) ==  1.3271599928087172e-24
end
