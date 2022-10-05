struct Nihalf
    n::Rational
    eps::Float64
end

@testset "lihalf" begin
    nis = (
        Nihalf(-20//2, 1e-09),
        Nihalf(-18//2, 1e-10),
        Nihalf(-16//2, 1e-10),
        Nihalf(-14//2, 1e-12),
        Nihalf(-12//2, 1e-12),
        Nihalf(-10//2, 1e-10),
        Nihalf( -8//2, 1e-13),
        Nihalf( -6//2, 1e-13),
        Nihalf( -4//2, 1e-13),
        Nihalf( -2//2, 1e-14),
        Nihalf(  0//2, 1e-14),
        Nihalf(  1//2, 1e-11), # TODO: increase
        Nihalf(  2//2, 1e-14),
        Nihalf(  3//2, 1e-11), # TODO: increase
        Nihalf(  4//2, 1e-14),
        Nihalf(  5//2, 1e-11), # TODO: increase
        Nihalf(  6//2, 1e-14),
        Nihalf(  8//2, 1e-14),
        Nihalf( 10//2, 1e-14),
        Nihalf( 12//2, 1e-14),
        # Nihalf( 21//2, 1e-10), # TODO: increase
        # Nihalf( 41//2, 1e-10), # TODO: increase
    )

    for ni in nis
        n = ni.n
        if denominator(n) == 1
            k = 2*numerator(n)
        elseif denominator(n) == 2
            k = numerator(n)
        else
            throw(DomainError(n, "only half-integer orders are allowed"))
        end

        if iseven(k)
            cmpl_data = read_from(joinpath(@__DIR__, "data", "Li$(k÷2).txt"))
            real_data = filter_real(cmpl_data)
            test_function_on_data(z -> PolyLog.reli(n, z), real_data, ni.eps, ni.eps)
            test_function_on_data(z -> PolyLog.li(n, z), cmpl_data, ni.eps, ni.eps)
        else
            cmpl_data = read_from(joinpath(@__DIR__, "data", "Li$(k)half.txt"))
            real_data = filter_real(cmpl_data)

            # remove Li_{1/2}(<≈1), which is Inf
            if n == 1//2
                real_data = hcat(
                    [real_data[i,1] for i in 1:size(real_data, 1) if !(real_data[i,1] <= one(Float64) && one(Float64) - real_data[i,1] < eps(Float64))],
                    [real_data[i,2] for i in 1:size(real_data, 1) if !(real_data[i,1] <= one(Float64) && one(Float64) - real_data[i,1] < eps(Float64))]
                )
                cmpl_data = hcat(
                    [cmpl_data[i,1] for i in 1:size(cmpl_data, 1) if !(abs(one(Float64) - cmpl_data[i,1]) < eps(Float64))],
                    [cmpl_data[i,2] for i in 1:size(cmpl_data, 1) if !(abs(one(Float64) - cmpl_data[i,1]) < eps(Float64))]
                )
            end

            test_function_on_data(z -> PolyLog.reli(n, z), real_data, ni.eps, ni.eps)
            test_function_on_data(z -> PolyLog.li(n, z), cmpl_data, ni.eps, ni.eps)
        end

        @test PolyLog.reli(n, 0.0) == zero(Float64)
        @test PolyLog.reli(n, 0.0f0) == zero(Float32)
        @test PolyLog.reli(n, Float16(0.0)) == zero(Float16)
        @test PolyLog.reli(n, 0//1) == 0
        @test PolyLog.reli(n, 0) == 0

        if n == 1//2
            @test PolyLog.reli(n, 1.0) == Inf
            @test PolyLog.reli(n, 1.0f0) == Inf
            @test PolyLog.reli(n, Float16(1.0)) == Inf
            @test PolyLog.reli(n, 1//1) == Inf
            @test PolyLog.reli(n, 1) == Inf
        else
            zeta = PolyLog.zetahalf(k)

            @test PolyLog.reli(n, 1.0) == zeta
            @test PolyLog.reli(n, 1.0f0) ≈ zeta
            @test PolyLog.reli(n, Float16(1.0)) ≈ zeta
            @test PolyLog.reli(n, 1//1) ≈ zeta
            @test PolyLog.reli(n, 1) ≈ zeta
        end

        @test isnan(PolyLog.reli(n, NaN))
        @test isinf(PolyLog.reli(n, Inf))
    end
end
