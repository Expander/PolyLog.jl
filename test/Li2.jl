@testset "li2" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li2.txt"), BigFloat)
    real_data = filter_real(cmpl_data)

    setprecision(BigFloat, MAX_BINARY_DIGITS) do
        ep = 10*eps(BigFloat)
        # @todo(alex): test imaginary part, too
        test_function_on_data(z -> real(PolyLog.li2(z)), real_data, ep, ep)
        test_function_on_data(PolyLog.reli2, real_data, ep, ep)
    end

    test_function_on_data(PolyLog.li2  , map(ComplexF64, cmpl_data), 1e-14, 1e-14)
    test_function_on_data(PolyLog.reli2, map(Float64   , real_data), 1e-14, 1e-14)

    test_function_on_data(PolyLog.li2  , map(ComplexF32, cmpl_data), 1e-6, 1e-6)
    test_function_on_data(PolyLog.reli2, map(Float32   , real_data), 1e-6, 1e-6)

    test_function_on_data(PolyLog.li2  , filter_ComplexF16(map(ComplexF16, cmpl_data)), 1e-2, 1e-2)
    test_function_on_data(PolyLog.reli2, map(Float16, real_data), 1e-2, 1e-2)

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

    # test value that causes overflow if squared
    @test real(PolyLog.li2(1e300 + 1im)) ≈ real(-238582.12510339421 - 2170.13532372464im) rtol=eps(Float64)
end
