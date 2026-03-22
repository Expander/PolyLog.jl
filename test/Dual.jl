if isdefined(Base, :get_extension)
    import ForwardDiff

    function li_dual(n::Integer, z::Complex)
        z_dual = complex(ForwardDiff.Dual(real(z)), ForwardDiff.Dual(imag(z)))
        PolyLog.li(n, z_dual)
    end


    @testset "Dual" begin
        z = 0.5 + 0.8im
        @test li_dual(2, z) == PolyLog.li(2, z)
    end
end
