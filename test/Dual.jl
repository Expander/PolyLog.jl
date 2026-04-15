if isdefined(Base, :get_extension)
    import ForwardDiff

    to_dual(z::Real) = ForwardDiff.Dual(z)
    to_dual(z::Complex) = Complex(ForwardDiff.Dual(real(z)), ForwardDiff.Dual(imag(z)))

    @testset "Dual" begin
        z = 0.5 + 0.8im

        @test PolyLog.reli1(to_dual(real(z))) ≈ PolyLog.reli1(real(z)) rtol=1e-14
        @test PolyLog.reli2(to_dual(real(z))) ≈ PolyLog.reli2(real(z)) rtol=1e-14
        @test PolyLog.reli3(to_dual(real(z))) ≈ PolyLog.reli3(real(z)) rtol=1e-14
        @test PolyLog.reli4(to_dual(real(z))) ≈ PolyLog.reli4(real(z)) rtol=1e-14

        @test PolyLog.li0(to_dual(z)) ≈ PolyLog.li0(z) rtol=1e-14
        @test PolyLog.li1(to_dual(z)) ≈ PolyLog.li1(z) rtol=1e-14
        @test PolyLog.li2(to_dual(z)) ≈ PolyLog.li2(z) rtol=1e-14

        for n in -10:10
            @test PolyLog.reli(n, to_dual(real(z))) ≈ PolyLog.reli(n, real(z)) rtol=1e-14
            @test PolyLog.li(n, to_dual(z)) ≈ PolyLog.li(n, z) rtol=1e-14
        end
    end
end
