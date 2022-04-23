@testset "li6" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li6.txt"))
    test_function_on_data(z -> PolyLog.li6(z), cmpl_data, 1e-14, 1e-14)
end
