@testset "li2" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li2.txt"))
    real_data = filter_real(cmpl_data)

    test_function_on_data(z -> PolyLog.li2(z)  , cmpl_data, 1e-14, 1e-14)
    test_function_on_data(z -> PolyLog.reli2(z), real_data, 1e-14, 1e-14)

    test_function_on_data(z -> PolyLog.li2(ComplexF32(z)), cmpl_data, 1e-6, 1e-6)
    test_function_on_data(z -> PolyLog.reli2(Float32(z) ), real_data, 1e-6, 1e-6)

    test_function_on_data(z -> PolyLog.li2(ComplexF16(z)), filter_ComplexF16(cmpl_data), 1e-2, 1e-2)
    test_function_on_data(z -> PolyLog.reli2(Float16(z) ), real_data, 1e-2, 1e-2)

    zeta2 = 1.6449340668482264

    @test PolyLog.reli2(1.0) == zeta2
    @test PolyLog.reli2(1.0f0) ≈ zeta2
    @test PolyLog.reli2(Float16(1.0)) ≈ zeta2
    @test PolyLog.reli2(1//1) == zeta2
    @test PolyLog.reli2(1) ≈ zeta2

    @test PolyLog.li2(1.0) == zeta2
    @test PolyLog.li2(1.0f0) ≈ zeta2
    @test PolyLog.li2(Float16(1.0)) ≈ zeta2
    @test PolyLog.li2(1//1) == zeta2
    @test PolyLog.li2(1) ≈ zeta2

    @test PolyLog.li2(1.0 + 0.0im) == zeta2
    @test PolyLog.li2(1.0f0 + 0.0f0im) ≈ zeta2
    @test PolyLog.li2(ComplexF16(1.0 + 0.0im)) ≈ zeta2
    @test PolyLog.li2(1//1 + 0//1im) ≈ zeta2
    @test PolyLog.li2(1 + 0im) ≈ zeta2

    @test real(PolyLog.li2(-1.08371e-08 + 1.32716e-24*im)) == -1.0837099970639316e-08
    @test imag(PolyLog.li2(-1.08371e-08 + 1.32716e-24*im)) ==  1.3271599928087172e-24
end
