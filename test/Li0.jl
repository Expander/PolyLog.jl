@testset "li0" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li0.txt"))
    real_data = filter_real(cmpl_data)

    test_function_on_data(z -> PolyLog.li0(z), cmpl_data, 1e-14, 1e-14)
    test_function_on_data(z -> PolyLog.li0(z), real_data, 1e-14, 1e-14)

    test_cmpl_function_on_data(z -> PolyLog.li0(ComplexF32(z)), joinpath(@__DIR__, "data", "Li0.txt"), 1e-6, 1e-6)
    test_real_function_on_data(z -> PolyLog.li0(Float32(z)   ), joinpath(@__DIR__, "data", "Li0.txt"), 1e-6, 1e-6)

    test_cmpl_function_on_data(z -> PolyLog.li0(ComplexF16(z)), joinpath(@__DIR__, "data", "Li0.txt"), 1e-2, 1e-2)
    test_real_function_on_data(z -> PolyLog.li0(Float16(z)   ), joinpath(@__DIR__, "data", "Li0.txt"), 1e-2, 1e-2)

    @test PolyLog.li0(1.0) == Inf
    @test isnan(PolyLog.li0(1.0 + 0.0im))
end
