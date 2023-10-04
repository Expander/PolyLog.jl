@testset "li0" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li0.txt"), BigFloat)
    real_data = filter_real(cmpl_data)

    setprecision(BigFloat, MAX_BINARY_DIGITS) do
        ep = 10*eps(BigFloat)
        test_function_on_data(z -> PolyLog.li0(z), cmpl_data, ep, ep)
    end

    test_function_on_data(PolyLog.li0, map(ComplexF64, cmpl_data), 1e-14, 1e-14)
    test_function_on_data(PolyLog.li0, map(Float64   , real_data), 1e-14, 1e-14)

    test_function_on_data(PolyLog.li0, map(ComplexF32, cmpl_data), 1e-6, 1e-6)
    test_function_on_data(PolyLog.li0, map(Float32   , real_data), 1e-6, 1e-6)

    test_function_on_data(PolyLog.li0, map(ComplexF16, cmpl_data), 1e-2, 1e-2)
    test_function_on_data(PolyLog.li0, map(Float16   , real_data), 1e-2, 1e-2)

    @test PolyLog.li0(2.0) ≈ -2.0
    @test PolyLog.li0(2.0f0) ≈ -2.0f0
    @test PolyLog.li0(Float16(2.0)) ≈ Float16(-2.0)
    @test PolyLog.li0(2//1) == -2//1
    @test PolyLog.li0(2) ≈ -2.0

    @test PolyLog.li0(2.0 + 0.0im) == -2.0
    @test PolyLog.li0(2.0f0 + 0.0f0im) ≈ -2.0f0
    @test PolyLog.li0(ComplexF16(2.0 + 0.0im)) ≈ Float16(-2.0)
    @test PolyLog.li0(2//1 + 0//1im) ≈ -2.0
    @test PolyLog.li0(2 + 0im) ≈ -2.0

    @test PolyLog.li0(1.0) == Inf
    @test isnan(PolyLog.li0(1.0 + 0.0im))
end
