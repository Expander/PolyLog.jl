@testset "eta" begin
    @test typeof(PolyLog.neg_eta(1, Float16)) == Float16
    @test typeof(PolyLog.neg_eta(1, Float32)) == Float32
    @test typeof(PolyLog.neg_eta(1, Float64)) == Float64
    @test typeof(PolyLog.neg_eta(1, BigFloat)) == BigFloat

    @test PolyLog.neg_eta(-222, Float64) ==  0.0
    @test PolyLog.neg_eta(-221, Float64) == -Inf
    @test PolyLog.neg_eta(-220, Float64) ==  0.0
    @test PolyLog.neg_eta(-219, Float64) == Inf
    @test PolyLog.neg_eta(-218, Float64) ==  0.0
    @test PolyLog.neg_eta(-217, Float64) == -1.8184610414701105e306
    @test PolyLog.neg_eta(  -4, Float64) ==  0.0
    @test PolyLog.neg_eta(  -3, Float64) ==  1/8
    @test PolyLog.neg_eta(  -2, Float64) ==  0.0
    @test PolyLog.neg_eta(  -1, Float64) == -0.25
    @test PolyLog.neg_eta(   0, Float64) == -0.5
    @test PolyLog.neg_eta(   1, Float64) == -0.69314718055994531
    @test PolyLog.neg_eta(   2, Float64) == -0.82246703342411322
    @test PolyLog.neg_eta(  52, Float64) == -0.9999999999999998
    @test PolyLog.neg_eta(  53, Float64) == -0.9999999999999999
    @test PolyLog.neg_eta(  54, Float64) == -0.9999999999999999
    @test PolyLog.neg_eta(  55, Float64) == -1.0
    @test PolyLog.neg_eta(  56, Float64) == -1.0

    @test PolyLog.neg_eta(-2, BigFloat) == BigFloat("0")
    @test PolyLog.neg_eta(-1, BigFloat) == (big(2)^2 - one(BigFloat))*PolyLog.zeta(-1, BigFloat)
    @test PolyLog.neg_eta( 0, BigFloat) == BigFloat("-0.5")
    @test PolyLog.neg_eta( 1, BigFloat) == -log(big(2))
    @test PolyLog.neg_eta( 2, BigFloat) == (inv(big(2)) - one(BigFloat))*PolyLog.zeta(2, BigFloat)
end
