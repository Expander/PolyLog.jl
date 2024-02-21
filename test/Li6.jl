@testset "li6" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li6.txt"), BigFloat)

    test_function_on_data(PolyLog.li6, map(ComplexF64, cmpl_data), 1e-14, 1e-14)
    test_function_on_data(PolyLog.li6, map(ComplexF32, cmpl_data), 1e-6, 1e-6)
    test_function_on_data(PolyLog.li6, map(ComplexF16, cmpl_data), 1e-2, 1e-2)

    # test signbit for 0.0 and -0.0 arguments
    for T in (Float16, Float32, Float64)
        @test  signbit(real(PolyLog.li6(T(-0.0))))
        @test !signbit(imag(PolyLog.li6(T(-0.0))))
        @test !signbit(real(PolyLog.li6(T( 0.0))))
        @test !signbit(imag(PolyLog.li6(T( 0.0))))
        @test !signbit(real(PolyLog.li6(Complex{T}(0.0, 0.0))))
        @test !signbit(imag(PolyLog.li6(Complex{T}(0.0, 0.0))))
        @test  signbit(real(PolyLog.li6(Complex{T}(-0.0, 0.0))))
        @test !signbit(imag(PolyLog.li6(Complex{T}(-0.0, 0.0))))
        @test !signbit(real(PolyLog.li6(Complex{T}(0.0, -0.0))))
        @test  signbit(imag(PolyLog.li6(Complex{T}(0.0, -0.0))))
        @test  signbit(real(PolyLog.li6(Complex{T}(-0.0, -0.0))))
        @test  signbit(imag(PolyLog.li6(Complex{T}(-0.0, -0.0))))
    end

    zeta6 = 1.0173430619844491

    @test PolyLog.li6(1.0 + 0.0im) == zeta6
    @test PolyLog.li6(1.0f0 + 0.0f0im) ≈ zeta6
    @test PolyLog.li6(ComplexF16(1.0 + 0.0im)) ≈ zeta6
    @test PolyLog.li6(1//1 + 0//1im) ≈ zeta6
    @test PolyLog.li6(1 + 0im) ≈ zeta6

    # test value that causes overflow if squared
    @test PolyLog.li6(1e300 + 1im) ≈ -1.5086876165613597e14 + 4.11768711823317e12im rtol=eps(Float64)
    @test PolyLog.li6(1.0 + 1e300im) ≈ -1.5090387516918862e14 + 2.0589500211678e12im rtol=eps(Float64)
end
