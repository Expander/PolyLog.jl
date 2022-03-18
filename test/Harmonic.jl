@testset "harmonic" begin
    @test_throws DomainError PolyLog.harmonic(-1)
    @test_throws DomainError PolyLog.harmonic(0)
    @test PolyLog.harmonic( 1) == 1.0
    @test PolyLog.harmonic( 2) == 3/2
    @test PolyLog.harmonic( 3) == 11/6
    @test PolyLog.harmonic( 4) ≈ 25/12              rtol=1e-14
    @test PolyLog.harmonic(19) ≈ 3.5477396571436819 rtol=1e-14
    @test PolyLog.harmonic(20) ≈ 3.5977396571436819 rtol=1e-14
    @test PolyLog.harmonic(21) ≈ 3.6453587047627295 rtol=1e-14
end
