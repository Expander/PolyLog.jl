@testset "li2" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li2.txt"), BigFloat)
    real_data = filter_real(cmpl_data)

    setprecision(BigFloat, MAX_BINARY_DIGITS) do
        ep = 30*eps(BigFloat)
        test_function_on_data(PolyLog.li2  , cmpl_data, ep, ep)
        ep = 10*eps(BigFloat)
        test_function_on_data(PolyLog.reli2, real_data, ep, ep)
    end

    test_function_on_data(PolyLog.li2  , map(ComplexF64, cmpl_data), 1e-14, 1e-14)
    test_function_on_data(PolyLog.reli2, map(Float64   , real_data), 1e-14, 1e-14)

    test_function_on_data(PolyLog.li2  , map(ComplexF32, cmpl_data), 1e-6, 1e-6)
    test_function_on_data(PolyLog.reli2, map(Float32   , real_data), 1e-6, 1e-6)

    test_function_on_data(PolyLog.li2  , filter_ComplexF16(map(ComplexF16, cmpl_data)), 1e-2, 1e-2)
    test_function_on_data(PolyLog.reli2, map(Float16, real_data), 1e-2, 1e-2)

    # test signbit for 0.0 and -0.0 arguments
    for T in (Float16, Float32, Float64, BigFloat)
        @test  signbit(PolyLog.reli2(T(-0.0)))
        @test !signbit(PolyLog.reli2(T( 0.0)))
        @test  signbit(real(PolyLog.li2(T(-0.0))))
        @test !signbit(imag(PolyLog.li2(T(-0.0))))
        @test !signbit(real(PolyLog.li2(T( 0.0))))
        @test !signbit(imag(PolyLog.li2(T( 0.0))))
        @test !signbit(real(PolyLog.li2(Complex{T}(0.0, 0.0))))
        @test !signbit(imag(PolyLog.li2(Complex{T}(0.0, 0.0))))
        @test  signbit(real(PolyLog.li2(Complex{T}(-0.0, 0.0))))
        @test !signbit(imag(PolyLog.li2(Complex{T}(-0.0, 0.0))))
        @test !signbit(real(PolyLog.li2(Complex{T}(0.0, -0.0))))
        @test  signbit(imag(PolyLog.li2(Complex{T}(0.0, -0.0))))
        @test  signbit(real(PolyLog.li2(Complex{T}(-0.0, -0.0))))
        @test  signbit(imag(PolyLog.li2(Complex{T}(-0.0, -0.0))))
    end

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

    # test value close to zero
    @test real(PolyLog.li2(4.831285545908206e-6 + 0.004396919500211628im)) ≈ -1.94166578202937687444628936853e-9 rtol=1e-12
    @test imag(PolyLog.li2(4.831285545908206e-6 + 0.004396919500211628im)) ≈  0.00439692067657240512726530759719387623 rtol=1e-14

    # test value that causes overflow if squared
    @test PolyLog.li2(1e300 + 1im) ≈ -238582.12510339421 + 2170.13532372464im rtol=eps(Float64)
    @test PolyLog.li2(1.0 + 1e300im) ≈ -238585.82620504462 + 1085.06766186232im rtol=eps(Float64)

    #ForwardDiff Test
    @test ForwardDiff.derivative(PolyLog.reli2,float(pi)) == PolyLog.reli1(float(pi))/float(pi)
end
