@testset "li0" begin
    test_cmpl_function_on_data(z -> PolyLog.li0(z), joinpath(@__DIR__, "data", "Li0.txt"), 1e-14, 1e-14)
    test_real_function_on_data(z -> PolyLog.li0(z), joinpath(@__DIR__, "data", "Li0.txt"), 1e-14, 1e-14)

    @test PolyLog.li0(1.0) == Inf
    @test isnan(PolyLog.li0(1.0 + 0.0im))
end
