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

        c64_data = read_from(joinpath(@__DIR__, "data", "Li$(n).txt"))
        f64_data = filter_real(c64_data)

        # test reli(n, z)
        for T in (Float16, Float32, Float64)
            ep = ni.eps*eps(T)/eps(Float64)
            test_function_on_data(z -> PolyLog.reli(n, T(z)), map(T, f64_data), ep, ep)
        end

        # test li(n, z)
        for T in (Float32, Float64)
            ep = ni.eps*eps(T)/eps(Float64)
            test_function_on_data(z -> PolyLog.li(n, z), map(Complex{T}, c64_data), ep, ep)
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

    setprecision(BigFloat, max_BigFloat_decimal_digits) do
        @test PolyLog.reli(5, BigFloat("0.5")) ≈ BigFloat("0.5084005792422687074591088492585899413195") atol=1e-39 rtol=1e-39
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
