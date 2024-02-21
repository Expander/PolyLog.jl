struct Ni
    n::Integer
    eps::Float64
end

@testset "li" begin
    nis = (
        Ni(    -10, 1e-09),
        Ni(     -9, 1e-10),
        Ni(     -8, 1e-10),
        Ni(     -7, 1e-12),
        Ni(     -6, 1e-12),
        Ni(     -5, 1e-10),
        Ni(     -4, 1e-13),
        Ni(     -3, 1e-13),
        Ni(     -2, 1e-13),
        Ni(     -1, 1e-14),
        Ni(      0, 1e-14),
        Ni(      1, 1e-14),
        Ni(      2, 1e-14),
        Ni(      3, 1e-14),
        Ni(      4, 1e-14),
        Ni(      5, 1e-14),
        Ni(      6, 1e-14),
        Ni(    100, 1e-14),
        Ni(1000000, 1e-14)
    )

    for ni in nis
        n = ni.n

        complex_data = read_from(joinpath(@__DIR__, "data", "Li$(n).txt"), BigFloat)
        real_data = filter_real(complex_data)

        # test reli(n, z)
        setprecision(BigFloat, MAX_BINARY_DIGITS) do
            for T in (Float16, Float32, Float64, BigFloat)
                ep = ni.eps*eps(T)/eps(Float64)
                for TN in (Int8, Int16, Int32, Int64, Int128)
                    (n > typemax(TN) || n < typemin(TN)) && continue
                    test_function_on_data(z -> PolyLog.reli(TN(n), z), map(T, real_data), ep, ep)
                end
            end
        end

        # test li(n, z)
        setprecision(BigFloat, MAX_BINARY_DIGITS) do
            for T in (Float32, Float64, BigFloat)
                T == BigFloat && (n < 0 || n > 2) && continue # tests take too long
                ep = ni.eps*eps(T)/eps(Float64)
                for TN in (Int8, Int16, Int32, Int64, Int128)
                    (n > typemax(TN) || n < typemin(TN)) && continue
                    test_function_on_data(z -> PolyLog.li(TN(n), z), map(Complex{T}, complex_data), ep, ep)
                end
            end
        end

        # test signbit for 0.0 and -0.0 arguments
        for T in (Float16, Float32, Float64)
            @test  signbit(real(PolyLog.li(n, T(-0.0))))
            @test !signbit(imag(PolyLog.li(n, T(-0.0))))
            @test !signbit(real(PolyLog.li(n, T( 0.0))))
            @test !signbit(imag(PolyLog.li(n, T( 0.0))))
        end

        zeta = PolyLog.zeta(n, Float64)

        @test PolyLog.reli(n, 1.0) == zeta
        @test PolyLog.reli(n, 1.0f0) ≈ zeta
        @test PolyLog.reli(n, Float16(1.0)) ≈ zeta
        @test PolyLog.reli(n, 1//1) ≈ zeta
        @test PolyLog.reli(n, 1) ≈ zeta
        @test PolyLog.reli(n, BigFloat("1.0")) == PolyLog.zeta(n, BigFloat)

        @test PolyLog.li(n, 1.0) == zeta
        @test PolyLog.li(n, 1.0f0) ≈ zeta
        @test PolyLog.li(n, Float16(1.0)) ≈ zeta
        @test PolyLog.li(n, 1//1) ≈ zeta
        @test PolyLog.li(n, 1) ≈ zeta
        @test PolyLog.li(n, BigFloat("1.0")) == PolyLog.zeta(n, BigFloat)

        @test PolyLog.li(n, 1.0 + 0.0im) == zeta
        @test PolyLog.li(n, 1.0f0 + 0.0f0im) ≈ zeta
        @test PolyLog.li(n, ComplexF16(1.0 + 0.0im)) ≈ zeta
        @test PolyLog.li(n, 1//1 + 0//1im) ≈ zeta
        @test PolyLog.li(n, 1 + 0im) ≈ zeta
        @test PolyLog.li(n, BigFloat("1.0") + 0im) == PolyLog.zeta(n, BigFloat)
    end

    # value close to boundary between series 1 and 2 in arXiv:2010.09860
    @test PolyLog.li(-2, -0.50001)               ≈ -0.074072592582716422 atol=1e-14
    @test PolyLog.reli(-2, -0.50001)             ≈ -0.074072592582716422 atol=1e-14
    @test PolyLog.li(-2, BigFloat("-0.50001"))   ≈ -0.074072592582716422 atol=1e-14
    @test PolyLog.reli(-2, BigFloat("-0.50001")) ≈ -0.074072592582716422 atol=1e-14

    # value sensitive to proper treatment of 0.0 vs -0.0 in imag(z)
    z = 1.5 + 0.0im
    @test PolyLog.li(10,  z) ≈ 1.5022603281703005298 - 2.56429642116111388671e-9im atol=1e-14 rtol=1e-14
    @test PolyLog.li(10, -z) ≈ -1.4978556954869267594 atol=1e-14

    # test value that causes overflow if squared
    @test PolyLog.li(7, 1e300 + 1im) ≈ -1.4886831990993457e16 + 4.74066248802866e14im rtol=2*eps(Float64)
    @test PolyLog.li(7, 1.0 + 1e300im) ≈ -1.489168315226607e16 + 2.3705150998401e14im rtol=2*eps(Float64)

    @test isnan(PolyLog.reli(10, NaN))
    @test isnan(PolyLog.reli(10, BigFloat(NaN)))
    @test isinf(PolyLog.reli(10, Inf))
    @test isinf(PolyLog.reli(10, BigFloat(Inf)))
    @test isnan(PolyLog.li(10, NaN))
    @test isnan(PolyLog.li(10, BigFloat(NaN)))
    @test isinf(PolyLog.li(10, Inf))
    @test isinf(PolyLog.li(10, BigFloat(Inf)))

    # test type stability
    for T in (Float16, Float32, Float64, BigFloat)
        for n in -10:10
            for x in (NaN, Inf, -2, -1, 0, 1, 2, 3)
                @test typeof(PolyLog.li(n, convert(Complex{T}, x))) == Complex{T}
                @test typeof(PolyLog.reli(n, convert(T, x))) == T
            end
        end
    end
end
