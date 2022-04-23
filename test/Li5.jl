@testset "li5" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li5.txt"))
    test_function_on_data(z -> PolyLog.li5(z), cmpl_data, 1e-14, 1e-14)
end
