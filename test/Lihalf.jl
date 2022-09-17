struct Nihalf
    n::Integer
    eps::Float64
end

@testset "lihalf" begin
    nis = (
        Ni(-20, 1e-09),
        Ni(-18, 1e-10),
        Ni(-16, 1e-10),
        Ni(-14, 1e-12),
        Ni(-12, 1e-12),
        Ni(-10, 1e-10),
        Ni( -8, 1e-13),
        Ni( -6, 1e-13),
        Ni( -4, 1e-13),
        Ni( -2, 1e-14),
        Ni(  0, 1e-14),
        Ni(  1, 1e-14),
        Ni(  2, 1e-14),
        Ni(  3, 1e-14),
        Ni(  4, 1e-14),
        Ni(  5, 1e-14),
        Ni(  6, 1e-14),
        Ni(  8, 1e-14),
        Ni( 10, 1e-14),
        Ni( 12, 1e-14),
    )

    for ni in nis
        n    = ni.n
        eps  = ni.eps

        if iseven(n)
            cmpl_data = read_from(joinpath(@__DIR__, "data", "Li$(n÷2).txt"))
            real_data = filter_real(cmpl_data)
            test_function_on_data(z -> PolyLog.relihalf(n, z), real_data, eps, eps)
        else
            cmpl_data = read_from(joinpath(@__DIR__, "data", "Li$(n)half.txt"))
            real_data = filter_real(cmpl_data)

            @test PolyLog.relihalf(n, 0.0) == zero(Float64)
            @test PolyLog.relihalf(n, 0.0f0) == zero(Float32)
            @test PolyLog.relihalf(n, Float16(0.0)) == zero(Float16)
            @test PolyLog.relihalf(n, 0) == 0
        end

        # zeta = PolyLog.zeta(n)

        # @test PolyLog.relihalf(n, 1.0) == zeta
        # @test PolyLog.relihalf(n, 1.0f0) ≈ zeta
        # @test PolyLog.relihalf(n, Float16(1.0)) ≈ zeta
        # @test PolyLog.relihalf(n, 1//1) ≈ zeta
        # @test PolyLog.relihalf(n, 1) ≈ zeta

        @test isnan(PolyLog.relihalf(n, NaN))
        @test isinf(PolyLog.relihalf(n, Inf))
    end
end
