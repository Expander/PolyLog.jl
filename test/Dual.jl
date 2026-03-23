if isdefined(Base, :get_extension)
    import ForwardDiff

    to_dual(z::Complex) = Complex(ForwardDiff.Dual(real(z)), ForwardDiff.Dual(imag(z)))

    @testset "Dual" begin
        z = 0.5 + 0.8im

        for n in -10:10
            @test PolyLog.li(n, to_dual(z)) ≈ PolyLog.li(n, z) rtol=1e-14
        end
    end
end
