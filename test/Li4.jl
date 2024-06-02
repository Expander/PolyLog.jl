@testset "li4" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li4.txt"), BigFloat)
    real_data = filter_real(cmpl_data)

    test_function_on_data(PolyLog.li4  , map(ComplexF64, cmpl_data), 1e-14, 1e-14)
    test_function_on_data(PolyLog.reli4, map(Float64   , real_data), 1e-14, 1e-14)

    test_function_on_data(PolyLog.li4  , map(ComplexF32, cmpl_data), 1e-6, 1e-6)
    test_function_on_data(PolyLog.reli4, map(Float32   , real_data), 1e-6, 1e-6)

    test_function_on_data(PolyLog.li4  , map(ComplexF16, cmpl_data), 1e-2, 1e-2)
    test_function_on_data(PolyLog.reli4, map(Float16   , real_data), 1e-2, 1e-2)

    # test signbit for 0.0 and -0.0 arguments
    for T in (Float16, Float32, Float64)
        @test  signbit(PolyLog.reli4(T(-0.0)))
        @test !signbit(PolyLog.reli4(T( 0.0)))
        @test  signbit(real(PolyLog.li4(T(-0.0))))
        @test !signbit(imag(PolyLog.li4(T(-0.0))))
        @test !signbit(real(PolyLog.li4(T( 0.0))))
        @test !signbit(imag(PolyLog.li4(T( 0.0))))
        @test !signbit(real(PolyLog.li4(Complex{T}(0.0, 0.0))))
        @test !signbit(imag(PolyLog.li4(Complex{T}(0.0, 0.0))))
        @test  signbit(real(PolyLog.li4(Complex{T}(-0.0, 0.0))))
        @test !signbit(imag(PolyLog.li4(Complex{T}(-0.0, 0.0))))
        @test !signbit(real(PolyLog.li4(Complex{T}(0.0, -0.0))))
        @test  signbit(imag(PolyLog.li4(Complex{T}(0.0, -0.0))))
        @test  signbit(real(PolyLog.li4(Complex{T}(-0.0, -0.0))))
        @test  signbit(imag(PolyLog.li4(Complex{T}(-0.0, -0.0))))
    end

    zeta4 = 1.0823232337111382

    @test PolyLog.reli4(1.0) == zeta4
    @test PolyLog.reli4(1.0f0) ≈ zeta4
    @test PolyLog.reli4(Float16(1.0)) ≈ zeta4
    @test PolyLog.reli4(1//1) ≈ zeta4
    @test PolyLog.reli4(1) ≈ zeta4

    @test PolyLog.li4(1.0) == zeta4
    @test PolyLog.li4(1.0f0) ≈ zeta4
    @test PolyLog.li4(Float16(1.0)) ≈ zeta4
    @test PolyLog.li4(1//1) ≈ zeta4
    @test PolyLog.li4(1) ≈ zeta4

    @test PolyLog.li4(1.0 + 0.0im) == zeta4
    @test PolyLog.li4(1.0f0 + 0.0f0im) ≈ zeta4
    @test PolyLog.li4(ComplexF16(1.0 + 0.0im)) ≈ zeta4
    @test PolyLog.li4(1//1 + 0//1im) ≈ zeta4
    @test PolyLog.li4(1 + 0im) ≈ zeta4

    # test value that causes overflow if squared
    @test PolyLog.li4(1e300 + 1im) ≈ -9.4863817894708364e9 + 1.725875455850714e8im rtol=eps(Float64)
    @test PolyLog.li4(1.0 + 1e300im) ≈ -9.4872648206269765e9 + 8.62951114411071e7im rtol=eps(Float64)

    #ForwardDiff Test
    @test ForwardDiff.derivative(PolyLog.reli4,float(pi)) == reli3(float(pi))/float(pi)
end
