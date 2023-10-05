@testset "clog" begin
    @test PolyLog.clog(-1.0 + 0.0im) == log(-1.0 + 0.0im)
    @test PolyLog.clog( 1.0 + 0.0im) == log(1.0)
    @test PolyLog.clog( 2.0 + 0.0im) == log(2.0)
    @test PolyLog.clog( 1.0) == log(1.0)
    @test PolyLog.clog( 2.0) == log(2.0)

    # test case where imag(z) == -0.0
    z = 1.0 + 0.0im
    @test imag(PolyLog.clog(-z)) ≈  pi rtol=1e-15
    @test imag(log(-z))          ≈ -pi rtol=1e-15
end

@testset "ln_sqr" begin
    for x in -10:10
        for T in (Float16, Float32, Float64, BigFloat)
            @test PolyLog.ln_sqr(T(x)) == abs2(log(complex(T(x))))
        end
    end
end
