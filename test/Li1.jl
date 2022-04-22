@testset "li1" begin
    test_cmpl_function_on_data(z -> PolyLog.li1(z), joinpath(@__DIR__, "data", "Li1.txt"), 1e-14, 1e-14)
    test_real_function_on_data(z -> PolyLog.li1(z), joinpath(@__DIR__, "data", "Li1.txt"), 1e-14, 1e-14)

    @test PolyLog.li1(1.0) == Inf
    @test PolyLog.li1(1.0 + 0.0im) == Inf
end
