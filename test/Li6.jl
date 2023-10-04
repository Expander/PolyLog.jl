@testset "li6" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li6.txt"), BigFloat)

    test_function_on_data(z -> PolyLog.li6(ComplexF64(z)), cmpl_data, 1e-14, 1e-14)
    test_function_on_data(z -> PolyLog.li6(ComplexF32(z)), cmpl_data, 1e-6, 1e-6)
    test_function_on_data(z -> PolyLog.li6(ComplexF16(z)), cmpl_data, 1e-2, 1e-2)

    zeta6 = 1.0173430619844491

    @test PolyLog.li6(1.0 + 0.0im) == zeta6
    @test PolyLog.li6(1.0f0 + 0.0f0im) ≈ zeta6
    @test PolyLog.li6(ComplexF16(1.0 + 0.0im)) ≈ zeta6
    @test PolyLog.li6(1//1 + 0//1im) ≈ zeta6
    @test PolyLog.li6(1 + 0im) ≈ zeta6
end
