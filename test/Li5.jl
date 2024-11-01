@testset "li5" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li5.txt"), BigFloat)

    test_function_on_data(PolyLog.li5, map(ComplexF64, cmpl_data), 1e-14, 1e-14)
    test_function_on_data(PolyLog.li5, map(ComplexF32, cmpl_data), 1e-6, 1e-6)
    test_function_on_data(PolyLog.li5, map(ComplexF16, cmpl_data), 1e-2, 1e-2)

    # test signbit for 0.0 and -0.0 arguments
    for T in (Float16, Float32, Float64)
        @test  signbit(real(PolyLog.li5(T(-0.0))))
        @test !signbit(imag(PolyLog.li5(T(-0.0))))
        @test !signbit(real(PolyLog.li5(T( 0.0))))
        @test !signbit(imag(PolyLog.li5(T( 0.0))))
        @test !signbit(real(PolyLog.li5(Complex{T}(0.0, 0.0))))
        @test !signbit(imag(PolyLog.li5(Complex{T}(0.0, 0.0))))
        @test  signbit(real(PolyLog.li5(Complex{T}(-0.0, 0.0))))
        @test !signbit(imag(PolyLog.li5(Complex{T}(-0.0, 0.0))))
        @test !signbit(real(PolyLog.li5(Complex{T}(0.0, -0.0))))
        @test  signbit(imag(PolyLog.li5(Complex{T}(0.0, -0.0))))
        @test  signbit(real(PolyLog.li5(Complex{T}(-0.0, -0.0))))
        @test  signbit(imag(PolyLog.li5(Complex{T}(-0.0, -0.0))))
    end

    zeta5 = 1.0369277551433699

    @test PolyLog.li5(1.0 + 0.0im) == zeta5
    @test PolyLog.li5(1.0f0 + 0.0f0im) ≈ zeta5
    @test PolyLog.li5(ComplexF16(1.0 + 0.0im)) ≈ zeta5
    @test PolyLog.li5(1//1 + 0//1im) ≈ zeta5
    @test PolyLog.li5(1 + 0im) ≈ zeta5

    # test value that causes overflow if squared
    @test PolyLog.li5(1e300 + 1im) ≈ -1.3105197831948743e12 + 2.980481322754618e10im rtol=eps(Float64)
    @test PolyLog.li5(1.0 + 1e300im) ≈ -1.31072310968392418e12 + 1.490286896860219e10im rtol=eps(Float64)
end
