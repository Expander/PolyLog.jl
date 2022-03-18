@testset "digamma" begin
    @test PolyLog.digamma( 1) ≈ -0.57721566490153286 rtol=1e-14
    @test PolyLog.digamma( 2) ≈  0.42278433509846714 rtol=1e-14
    @test PolyLog.digamma( 3) ≈  0.92278433509846714 rtol=1e-14
    @test PolyLog.digamma(10) ≈  2.2517525890667211  rtol=1e-15
end
