@testset "harmonic" begin
    @test_throws DomainError PolyLog.harmonic(-1, Float64)
    @test_throws DomainError PolyLog.harmonic(-1, BigFloat)
    @test_throws DomainError PolyLog.harmonic(0, Float64)
    @test_throws DomainError PolyLog.harmonic(0, BigFloat)
    @test PolyLog.harmonic( 1, Float64) == 1.0
    @test PolyLog.harmonic( 2, Float64) == 3/2
    @test PolyLog.harmonic( 3, Float64) == 11/6
    @test PolyLog.harmonic( 4, Float64) ≈ 25/12              rtol=1e-14
    @test PolyLog.harmonic( 1, BigFloat) == BigFloat("1.0")
    @test PolyLog.harmonic( 2, BigFloat) == BigFloat("1.5")
    @test PolyLog.harmonic( 3, BigFloat) == BigFloat("11")/6
    @test PolyLog.harmonic( 4, BigFloat) ≈ BigFloat("25")/12 rtol=eps(BigFloat)
    @test PolyLog.harmonic(19, Float64) ≈ 3.5477396571436819 rtol=1e-14
    @test PolyLog.harmonic(20, Float64) ≈ 3.5977396571436819 rtol=1e-14
    @test PolyLog.harmonic(21, Float64) ≈ 3.6453587047627295 rtol=1e-14
end
