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
        data = open(readdlm, joinpath(@__DIR__, "data", "Li$(n).txt"))

        for r in 1:size(data, 1)
            row      = data[r, :]
            z        = row[1] + row[2]*1im
            expected = row[3] + row[4]*1im

            if imag(z) == 0.0
                @test PolyLog.li(n, real(z)) ≈ real(expected) atol=eps rtol=eps
            end
            if n <= 6
                @test PolyLog.li(n, z) ≈ expected atol=eps rtol=eps
            end
        end
    end

    # value close to boundary between series 1 and 2 in arXiv:2010.09860
    @test PolyLog.li(-2, -0.50001) ≈ -0.074072592582716422 atol=1e-14

    @test isnan(PolyLog.li(10, NaN))
    @test isinf(PolyLog.li(10, Inf))
end
