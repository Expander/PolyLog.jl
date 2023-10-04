@testset "li4" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li4.txt"), BigFloat)
    real_data = filter_real(cmpl_data)

    test_function_on_data(PolyLog.li4  , map(ComplexF64, cmpl_data), 1e-14, 1e-14)
    test_function_on_data(PolyLog.reli4, map(Float64   , real_data), 1e-14, 1e-14)

    test_function_on_data(PolyLog.li4  , map(ComplexF32, cmpl_data), 1e-6, 1e-6)
    test_function_on_data(PolyLog.reli4, map(Float32   , real_data), 1e-6, 1e-6)

    test_function_on_data(PolyLog.li4  , map(ComplexF16, cmpl_data), 1e-2, 1e-2)
    test_function_on_data(PolyLog.reli4, map(Float16   , real_data), 1e-2, 1e-2)

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
end
