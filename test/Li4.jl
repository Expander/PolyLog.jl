@testset "li4" begin
    cmpl_data = read_from(joinpath(@__DIR__, "data", "Li4.txt"))
    real_data = filter_real(cmpl_data)

    test_function_on_data(z -> PolyLog.li4(z), cmpl_data, 1e-14, 1e-14)
    test_function_on_data(z -> PolyLog.li4(z), real_data, 1e-14, 1e-14)
end
