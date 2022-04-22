@testset "li5" begin
    test_cmpl_function_on_data(z -> PolyLog.li5(z), joinpath(@__DIR__, "data", "Li5.txt"), 1e-14, 1e-14)
end
