struct Nihalf
    n::Integer
    eps::Float64
end

@testset "lihalf" begin
    nis = (
        Nihalf(-20, 1e-09),
        Nihalf(-18, 1e-10),
        Nihalf(-16, 1e-10),
        Nihalf(-14, 1e-12),
        Nihalf(-12, 1e-12),
        Nihalf(-10, 1e-10),
        Nihalf( -8, 1e-13),
        Nihalf( -6, 1e-13),
        Nihalf( -4, 1e-13),
        Nihalf( -2, 1e-14),
        Nihalf(  0, 1e-14),
        Nihalf(  1, 1e-11), # TODO: increase
        Nihalf(  2, 1e-14),
        Nihalf(  3, 1e-11), # TODO: increase
        Nihalf(  4, 1e-14),
        Nihalf(  5, 1e-11), # TODO: increase
        Nihalf(  6, 1e-14),
        Nihalf(  8, 1e-14),
        Nihalf( 10, 1e-14),
        Nihalf( 12, 1e-14),
        # Nihalf( 21, 1e-10), # TODO: increase
        # Nihalf( 41, 1e-10), # TODO: increase
    )

    for ni in nis
        n = ni.n

        if iseven(n)
            cmpl_data = read_from(joinpath(@__DIR__, "data", "Li$(n÷2).txt"))
            real_data = filter_real(cmpl_data)
            test_function_on_data(z -> PolyLog.relihalf(n, z), real_data, ni.eps, ni.eps)
            test_function_on_data(z -> PolyLog.lihalf(n, z), cmpl_data, ni.eps, ni.eps)
        else
            cmpl_data = read_from(joinpath(@__DIR__, "data", "Li$(n)half.txt"))
            real_data = filter_real(cmpl_data)

            # remove Li_{1/2}(<≈1), which is Inf
            if n == 1
                real_data = hcat(
                    [real_data[i,1] for i in 1:size(real_data, 1) if !(real_data[i,1] <= one(Float64) && one(Float64) - real_data[i,1] < eps(Float64))],
                    [real_data[i,2] for i in 1:size(real_data, 1) if !(real_data[i,1] <= one(Float64) && one(Float64) - real_data[i,1] < eps(Float64))]
                )
                cmpl_data = hcat(
                    [cmpl_data[i,1] for i in 1:size(cmpl_data, 1) if !(abs(one(Float64) - cmpl_data[i,1]) < eps(Float64))],
                    [cmpl_data[i,2] for i in 1:size(cmpl_data, 1) if !(abs(one(Float64) - cmpl_data[i,1]) < eps(Float64))]
                )
            end

            # test only values for which there is an implementation yet (TODO)
            # real_data = hcat(
            #     [real_data[i,1] for i in 1:size(real_data, 1) if abs(real_data[i,1]) < 0.75 || abs(log(Complex(real_data[i,1]))) < 2*pi],
            #     [real_data[i,2] for i in 1:size(real_data, 1) if abs(real_data[i,1]) < 0.75 || abs(log(Complex(real_data[i,1]))) < 2*pi]
            # )

            test_function_on_data(z -> PolyLog.relihalf(n, z), real_data, ni.eps, ni.eps)
            test_function_on_data(z -> PolyLog.lihalf(n, z), cmpl_data, ni.eps, ni.eps)
        end

        @test PolyLog.relihalf(n, 0.0) == zero(Float64)
        @test PolyLog.relihalf(n, 0.0f0) == zero(Float32)
        @test PolyLog.relihalf(n, Float16(0.0)) == zero(Float16)
        @test PolyLog.relihalf(n, 0//1) == 0
        @test PolyLog.relihalf(n, 0) == 0

        if n == 1
            @test PolyLog.relihalf(n, 1.0) == Inf
            @test PolyLog.relihalf(n, 1.0f0) == Inf
            @test PolyLog.relihalf(n, Float16(1.0)) == Inf
            @test PolyLog.relihalf(n, 1//1) == Inf
            @test PolyLog.relihalf(n, 1) == Inf
        else
            zeta = PolyLog.zetahalf(n)

            @test PolyLog.relihalf(n, 1.0) == zeta
            @test PolyLog.relihalf(n, 1.0f0) ≈ zeta
            @test PolyLog.relihalf(n, Float16(1.0)) ≈ zeta
            @test PolyLog.relihalf(n, 1//1) ≈ zeta
            @test PolyLog.relihalf(n, 1) ≈ zeta
        end

        @test isnan(PolyLog.relihalf(n, NaN))
        @test isinf(PolyLog.relihalf(n, Inf))
    end
end