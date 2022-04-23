@testset "li1" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li1.txt"))
    real_data = filter_real(cmpl_data)

    test_function_on_data(z -> PolyLog.li1(z), cmpl_data, 1e-14, 1e-14)
    test_function_on_data(z -> PolyLog.li1(z), real_data, 1e-14, 1e-14)

    test_function_on_data(z -> PolyLog.li1(ComplexF32(z)), cmpl_data, 1e-6, 1e-6)
    test_function_on_data(z -> PolyLog.li1(Float32(z)   ), real_data, 1e-6, 1e-6)

    test_function_on_data(z -> PolyLog.li1(ComplexF16(z)), cmpl_data, 1e-2, 1e-2)
    test_function_on_data(z -> PolyLog.li1(Float16(z)   ), real_data, 1e-2, 1e-2)

    @test PolyLog.li1(-1.0f0) ≈ -log(2.0f0)
    @test PolyLog.li1(-1//1) == -log(2.0)
    @test PolyLog.li1(-1) ≈ -log(2.0)
    @test PolyLog.li1(-1.0 + 0.0im) == -log(2.0)
    @test PolyLog.li1(-1.0f0 + 0.0f0im) ≈ -log(2.0f0)
    @test PolyLog.li1(-1//1 + 0//1im) ≈ -log(2.0)
    @test PolyLog.li1(-1 + 0im) ≈ -log(2.0)

    @test PolyLog.li1(1.0) == Inf
    @test PolyLog.li1(1.0 + 0.0im) == Inf
end
