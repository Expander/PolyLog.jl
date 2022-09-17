struct Nihalf
    n::Integer
    eps::Float64
end

@testset "lihalf" begin
    nis = (
        Ni(  0, 1e-14),
        Ni(  1, 1e-14),
        Ni(  3, 1e-14),
        Ni(  5, 1e-14),
    )

    for ni in nis
        n    = ni.n
        eps  = ni.eps

        cmpl_data = read_from(joinpath(@__DIR__, "data", "Li$(n)half.txt"))
        real_data = filter_real(cmpl_data)

        if n == 0
            test_function_on_data(z -> PolyLog.relihalf(n, z), real_data, eps, eps)
        end

        @test PolyLog.relihalf(n, 0.0) == zero(Float64)
        @test PolyLog.relihalf(n, 0.0f0) == zero(Float32)
        @test PolyLog.relihalf(n, Float16(0.0)) == zero(Float16)
        @test PolyLog.relihalf(n, 0) == 0

        # zeta = PolyLog.zeta(n)

        # @test PolyLog.relihalf(n, 1.0) == zeta
        # @test PolyLog.relihalf(n, 1.0f0) ≈ zeta
        # @test PolyLog.relihalf(n, Float16(1.0)) ≈ zeta
        # @test PolyLog.relihalf(n, 1//1) ≈ zeta
        # @test PolyLog.relihalf(n, 1) ≈ zeta
    end

    @test isnan(PolyLog.relihalf(10, NaN))
    @test isinf(PolyLog.relihalf(10, Inf))
end
