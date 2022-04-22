@testset "li6" begin
    test_cmpl_function_on_data(z -> PolyLog.li6(z), joinpath(@__DIR__, "data", "Li6.txt"), 1e-14)
end
