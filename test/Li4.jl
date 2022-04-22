@testset "li4" begin
    test_cmpl_function_on_data(z -> PolyLog.li4(z), joinpath(@__DIR__, "data", "Li4.txt"), 1e-14)
    test_real_function_on_data(z -> PolyLog.li4(z), joinpath(@__DIR__, "data", "Li4.txt"), 1e-14)
end
