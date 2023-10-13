@testset "li1" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li1.txt"), BigFloat)
    real_data = filter_real(cmpl_data)

    setprecision(BigFloat, MAX_BINARY_DIGITS) do
        ep = 10*eps(BigFloat)
        test_function_on_data(PolyLog.li1  , cmpl_data, ep, ep)
        test_function_on_data(PolyLog.reli1, real_data, ep, ep)
    end

    test_function_on_data(PolyLog.li1  , map(ComplexF64, cmpl_data), 1e-14, 1e-14)
    test_function_on_data(PolyLog.reli1, map(Float64   , real_data), 1e-14, 1e-14)

    test_function_on_data(PolyLog.li1  , map(ComplexF32, cmpl_data), 1e-6, 1e-6)
    test_function_on_data(PolyLog.reli1, map(Float32   , real_data), 1e-6, 1e-6)

    test_function_on_data(PolyLog.li1  , map(ComplexF16, cmpl_data), 1e-2, 1e-2)
    test_function_on_data(PolyLog.reli1, map(Float16   , real_data), 1e-2, 1e-2)

    @test PolyLog.reli1(-1.0) ≈ -log(2.0)
    @test PolyLog.reli1(-1.0f0) ≈ -log(2.0f0)
    @test PolyLog.reli1(Float16(-1.0)) ≈ Float16(-log(2.0))
    @test PolyLog.reli1(-1//1) == -log(2.0)
    @test PolyLog.reli1(-1) ≈ -log(2.0)
    @test PolyLog.reli1(1.0) == Inf

    @test PolyLog.li1(-1.0) ≈ -log(2.0)
    @test PolyLog.li1(-1.0f0) ≈ -log(2.0f0)
    @test PolyLog.li1(Float16(-1.0)) ≈ Float16(-log(2.0))
    @test PolyLog.li1(-1//1) == -log(2.0)
    @test PolyLog.li1(-1) ≈ -log(2.0)
    @test PolyLog.li1(1.0) == Inf

    @test PolyLog.li1(-1.0 + 0.0im) == -log(2.0)
    @test PolyLog.li1(-1.0f0 + 0.0f0im) ≈ -log(2.0f0)
    @test PolyLog.li1(ComplexF16(-1.0 + 0.0im)) ≈ Float16(-log(2.0))
    @test PolyLog.li1(-1//1 + 0//1im) ≈ -log(2.0)
    @test PolyLog.li1(-1 + 0im) ≈ -log(2.0)
    @test PolyLog.li1(1.0 + 0.0im) == Inf

    # test value that causes overflow if squared
    @test real(PolyLog.li1(1e300 + 1im)) ≈ real(-690.77552789821371 - 3.14159265358979im) rtol=eps(Float64)
end
