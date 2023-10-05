@testset "li5" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li5.txt"), BigFloat)

    test_function_on_data(PolyLog.li5, map(ComplexF64, cmpl_data), 1e-14, 1e-14)
    test_function_on_data(PolyLog.li5, map(ComplexF32, cmpl_data), 1e-6, 1e-6)
    test_function_on_data(PolyLog.li5, map(ComplexF16, cmpl_data), 1e-2, 1e-2)

    zeta5 = 1.0369277551433699

    @test PolyLog.li5(1.0 + 0.0im) == zeta5
    @test PolyLog.li5(1.0f0 + 0.0f0im) ≈ zeta5
    @test PolyLog.li5(ComplexF16(1.0 + 0.0im)) ≈ zeta5
    @test PolyLog.li5(1//1 + 0//1im) ≈ zeta5
    @test PolyLog.li5(1 + 0im) ≈ zeta5
end
