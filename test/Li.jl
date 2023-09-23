struct Ni
    n::Integer
    eps::Float64
end

const max_BigFloat_decimal_digits = ceil(Integer, 40*log(10)/log(2))

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
        n = ni.n

        complex_data = read_from(joinpath(@__DIR__, "data", "Li$(n).txt"), BigFloat)
        real_data = filter_real(complex_data)

        # test reli(n, z)
        for T in (Float16, Float32, Float64)
            ep = ni.eps*eps(T)/eps(Float64)
            test_function_on_data(z -> PolyLog.reli(n, T(z)), map(T, real_data), ep, ep)
        end

        # test reli(n, z) with BigFloat precision
        # @todo(alex): refactor
        setprecision(BigFloat, max_BigFloat_decimal_digits) do
            ep = 1e-39
            for i in 1:size(real_data, 1)
                z = real_data[i,1]
                if n > 4 && 0 < z && z < 3/4
                    x = PolyLog.reli(n, z)
                    y = real_data[i,2]
                    @test x ≈ y atol=ep rtol=ep
                end
            end
        end

        # test li(n, z)
        for T in (Float32, Float64)
            ep = ni.eps*eps(T)/eps(Float64)
            test_function_on_data(z -> PolyLog.li(n, z), map(Complex{T}, complex_data), ep, ep)
        end

        zeta = PolyLog.zeta(n)

        @test PolyLog.reli(n, 1.0) == zeta
        @test PolyLog.reli(n, 1.0f0) ≈ zeta
        @test PolyLog.reli(n, Float16(1.0)) ≈ zeta
        @test PolyLog.reli(n, 1//1) ≈ zeta
        @test PolyLog.reli(n, 1) ≈ zeta

        @test PolyLog.li(n, 1.0) == zeta
        @test PolyLog.li(n, 1.0f0) ≈ zeta
        @test PolyLog.li(n, Float16(1.0)) ≈ zeta
        @test PolyLog.li(n, 1//1) ≈ zeta
        @test PolyLog.li(n, 1) ≈ zeta

        @test PolyLog.li(n, 1.0 + 0.0im) == zeta
        @test PolyLog.li(n, 1.0f0 + 0.0f0im) ≈ zeta
        @test PolyLog.li(n, ComplexF16(1.0 + 0.0im)) ≈ zeta
        @test PolyLog.li(n, 1//1 + 0//1im) ≈ zeta
        @test PolyLog.li(n, 1 + 0im) ≈ zeta
    end

    # value close to boundary between series 1 and 2 in arXiv:2010.09860
    @test PolyLog.li(-2, -0.50001)   ≈ -0.074072592582716422 atol=1e-14
    @test PolyLog.reli(-2, -0.50001) ≈ -0.074072592582716422 atol=1e-14

    # value sensitive to proper treatment of 0.0 vs -0.0 in imag(z)
    z = 1.5 + 0.0im
    @test PolyLog.li(10,  z) ≈ 1.5022603281703005298 - 2.56429642116111388671e-9im atol=1e-14 rtol=1e-14
    @test PolyLog.li(10, -z) ≈ -1.4978556954869267594 atol=1e-14

    @test isnan(PolyLog.reli(10, NaN))
    @test isinf(PolyLog.reli(10, Inf))
    @test isnan(PolyLog.li(10, NaN))
    @test isinf(PolyLog.li(10, Inf))
end
