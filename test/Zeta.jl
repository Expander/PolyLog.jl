@testset "zeta" begin
    @test PolyLog.zeta(-263, Float64) == Inf
    @test PolyLog.zeta(-262, Float64) == 0.0
    @test PolyLog.zeta(-261, Float64) == -Inf
    @test PolyLog.zeta(-259, Float64) == 8.7601563446229215e306
    @test PolyLog.zeta(-257, Float64) == -5.1754977470366798e303
    @test PolyLog.zeta(-255, Float64) == 3.1055517596048927e300
    @test PolyLog.zeta(  -4, Float64) == 0.0
    @test PolyLog.zeta(  -3, Float64) == 1/120
    @test PolyLog.zeta(  -2, Float64) == 0.0
    @test PolyLog.zeta(  -1, Float64) == -1/12
    @test PolyLog.zeta(   0, Float64) == -0.5
    @test PolyLog.zeta(   1, Float64) == Inf
    @test PolyLog.zeta(   2, Float64) == pi^2/6
    @test PolyLog.zeta(  32, Float64) == 1.0000000002328312
    @test PolyLog.zeta(  33, Float64) == 1.0000000001164155
    @test PolyLog.zeta(  34, Float64) == 1.0000000000582077
    @test PolyLog.zeta(  35, Float64) == 1.0000000000291039

    @test PolyLog.zeta(-262, BigFloat) == BigFloat("0.0")
    @test PolyLog.zeta(  -4, BigFloat) == BigFloat("0.0")
    @test PolyLog.zeta(  -3, BigFloat) == big(1)/120
    @test PolyLog.zeta(  -2, BigFloat) == BigFloat("0.0")
    @test PolyLog.zeta(  -1, BigFloat) == -big(1)/12
    @test PolyLog.zeta(   0, BigFloat) == BigFloat("-0.5")
    @test PolyLog.zeta(   1, BigFloat) == BigFloat(Inf)
    @test PolyLog.zeta(   2, BigFloat) ≈ big(pi)^2/6       rtol=eps(BigFloat)
    @test PolyLog.zeta(   4, BigFloat) ≈ big(pi)^4/90      rtol=eps(BigFloat)
    @test PolyLog.zeta(   6, BigFloat) ≈ big(pi)^6/945     rtol=10*eps(BigFloat)
    @test PolyLog.zeta(   8, BigFloat) ≈ big(pi)^8/9450    rtol=10*eps(BigFloat)
    @test PolyLog.zeta(  10, BigFloat) ≈ big(pi)^10/93555  rtol=10*eps(BigFloat)
end
