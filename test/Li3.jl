@testset "li3" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li3.txt"), BigFloat)
    real_data = filter_real(cmpl_data)

    test_function_on_data(PolyLog.li3  , map(ComplexF64, cmpl_data), 1e-14, 1e-14)
    test_function_on_data(PolyLog.reli3, map(Float64   , real_data), 1e-14, 1e-14)

    test_function_on_data(PolyLog.li3  , map(ComplexF32, cmpl_data), 1e-6, 1e-6)
    test_function_on_data(PolyLog.reli3, map(Float32   , real_data), 1e-6, 1e-6)

    test_function_on_data(PolyLog.li3  , filter_ComplexF16(map(ComplexF16, cmpl_data)), 1e-2, 1e-2)
    test_function_on_data(PolyLog.reli3, map(Float16, real_data), 1e-2, 1e-2)

    # test signbit for 0.0 and -0.0 arguments
    for T in (Float16, Float32, Float64)
        @test  signbit(PolyLog.reli3(T(-0.0)))
        @test !signbit(PolyLog.reli3(T( 0.0)))
        @test  signbit(real(PolyLog.li3(T(-0.0))))
        @test !signbit(imag(PolyLog.li3(T(-0.0))))
        @test !signbit(real(PolyLog.li3(T( 0.0))))
        @test !signbit(imag(PolyLog.li3(T( 0.0))))
        @test !signbit(real(PolyLog.li3(Complex{T}(0.0, 0.0))))
        @test !signbit(imag(PolyLog.li3(Complex{T}(0.0, 0.0))))
        @test  signbit(real(PolyLog.li3(Complex{T}(-0.0, 0.0))))
        @test !signbit(imag(PolyLog.li3(Complex{T}(-0.0, 0.0))))
        @test !signbit(real(PolyLog.li3(Complex{T}(0.0, -0.0))))
        @test  signbit(imag(PolyLog.li3(Complex{T}(0.0, -0.0))))
        @test  signbit(real(PolyLog.li3(Complex{T}(-0.0, -0.0))))
        @test  signbit(imag(PolyLog.li3(Complex{T}(-0.0, -0.0))))
    end

    zeta3 = 1.2020569031595943

    @test PolyLog.reli3(1.0) == zeta3
    @test PolyLog.reli3(1.0f0) ≈ zeta3
    @test PolyLog.reli3(Float16(1.0)) ≈ zeta3
    @test PolyLog.reli3(1//1) ≈ zeta3
    @test PolyLog.reli3(1) ≈ zeta3

    @test PolyLog.reli3(0.0) == 0.0
    @test PolyLog.reli3(0.5) == 0.53721319360804020 # (-2*Pi^2*log(2) + 4*log(2)^3 + 21*zeta(3))/24
    @test PolyLog.reli3(-1.0) == -0.90154267736969571 # -3/4*zeta(3)

    @test PolyLog.li3(1.0) == zeta3
    @test PolyLog.li3(1.0f0) ≈ zeta3
    @test PolyLog.li3(Float16(1.0)) ≈ zeta3
    @test PolyLog.li3(1//1) ≈ zeta3
    @test PolyLog.li3(1) ≈ zeta3

    @test PolyLog.li3(1.0 + 0.0im) == zeta3
    @test PolyLog.li3(1.0f0 + 0.0f0im) ≈ zeta3
    @test PolyLog.li3(ComplexF16(1.0 + 0.0im)) ≈ zeta3
    @test PolyLog.li3(1//1 + 0//1im) ≈ zeta3
    @test PolyLog.li3(1 + 0im) ≈ zeta3

    @test PolyLog.li3(0.0) == 0.0
    @test PolyLog.li3(0.5) == 0.53721319360804020 # (-2*Pi^2*log(2) + 4*log(2)^3 + 21*zeta(3))/24
    @test PolyLog.li3(-1.0) == -0.90154267736969571 # -3/4*zeta(3)

    # test value that causes overflow if squared
    @test PolyLog.li3(1e300 + 1im) ≈ -5.4934049431527088e7 + 749538.186928224im rtol=eps(Float64)
    @test PolyLog.li3(1.0 + 1e300im) ≈ -5.4936606061973454e7 + 374771.031356405im rtol=eps(Float64)

    #ForwardDiff Test
    @test ForwardDiff.derivative(PolyLog.reli3,float(pi)) == reli2(float(pi))/float(pi)
end
