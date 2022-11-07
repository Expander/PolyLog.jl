@testset "digamma" begin
    @test_throws DomainError PolyLog.digamma(-1)
    @test_throws DomainError PolyLog.digamma(0)
    @test PolyLog.digamma(      1) ≈ -0.57721566490153286 rtol=1e-14
    @test PolyLog.digamma(      2) ≈  0.42278433509846714 rtol=1e-14
    @test PolyLog.digamma(      3) ≈  0.92278433509846714 rtol=1e-14
    @test PolyLog.digamma(      6) ≈  1.7061176684318005  rtol=1e-14
    @test PolyLog.digamma(      7) ≈  1.8727843350984671  rtol=1e-14
    @test PolyLog.digamma(     10) ≈  2.2517525890667211  rtol=1e-15
    @test PolyLog.digamma(    100) ≈  4.6001618527380874  rtol=1e-15
    @test PolyLog.digamma(   1000) ≈  6.9072551956488121  rtol=1e-15
    @test PolyLog.digamma(1000000) ≈  13.815510057964191  rtol=1e-15
end
