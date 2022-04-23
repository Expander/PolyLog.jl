struct Ni
    n::Integer
    eps::Float64
end

@testset "li" begin
    nis = (
        Ni(-10, 1e-09),
        Ni( -9, 1e-10),
        Ni( -8, 1e-10),
        Ni( -7, 1e-12),
        Ni( -6, 1e-12),
        Ni( -5, 1e-10),
        Ni( -4, 1e-13),
        Ni( -3, 1e-13),
        Ni( -2, 1e-13),
        Ni( -1, 1e-14),
        Ni(  0, 1e-14),
        Ni(  1, 1e-14),
        Ni(  2, 1e-14),
        Ni(  3, 1e-14),
        Ni(  4, 1e-14),
        Ni(  5, 1e-14),
        Ni(  6, 1e-14),
        Ni(100, 1e-14)
    )

    for ni in nis
        n    = ni.n
        eps  = ni.eps

        cmpl_data = read_from(joinpath(@__DIR__, "data", "Li$(n).txt"))
        real_data = filter_real(cmpl_data)

        test_function_on_data(z -> PolyLog.li(n, z), cmpl_data, eps, eps)
        test_function_on_data(z -> PolyLog.li(n, z), real_data, eps, eps)
    end

    # value close to boundary between series 1 and 2 in arXiv:2010.09860
    @test PolyLog.li(-2, -0.50001) ≈ -0.074072592582716422 atol=1e-14

    # value sensitive to proper treatment of 0.0 vs -0.0 in imag(z)
    z = 1.5 + 0.0im
    @test PolyLog.li(10,  z) ≈ 1.5022603281703005298 - 2.56429642116111388671e-9im atol=1e-14 rtol=1e-14
    @test PolyLog.li(10, -z) ≈ -1.4978556954869267594 atol=1e-14

    @test isnan(PolyLog.li(10, NaN))
    @test isinf(PolyLog.li(10, Inf))
end
