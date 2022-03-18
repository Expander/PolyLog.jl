@testset "eta" begin
    @test PolyLog.neg_eta(-222) ==  0.0
    @test PolyLog.neg_eta(-221) == -Inf
    @test PolyLog.neg_eta(-220) ==  0.0
    @test PolyLog.neg_eta(-219) == Inf
    @test PolyLog.neg_eta(-218) ==  0.0
    @test PolyLog.neg_eta(-217) == -1.8184610414701105e306
    @test PolyLog.neg_eta(  -4) ==  0.0
    @test PolyLog.neg_eta(  -3) ==  1/8
    @test PolyLog.neg_eta(  -2) ==  0.0
    @test PolyLog.neg_eta(  -1) == -0.25
    @test PolyLog.neg_eta(   1) == -0.69314718055994531
    @test PolyLog.neg_eta(   2) == -0.82246703342411322
    @test PolyLog.neg_eta(  52) == -0.9999999999999998
    @test PolyLog.neg_eta(  53) == -0.9999999999999999
    @test PolyLog.neg_eta(  54) == -0.9999999999999999
    @test PolyLog.neg_eta(  55) == -1.0
    @test PolyLog.neg_eta(  56) == -1.0
end
