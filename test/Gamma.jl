@testset "gamma" begin
    @test PolyLog.gammahalf(-341) == 0.0
    @test PolyLog.gammahalf(-339) == 5.6482208842233255e-306
    @test PolyLog.gammahalf(-3) == 2.3632718012073547
    @test PolyLog.gammahalf(-1) == -3.5449077018110321 # sqrt(pi)
    @test PolyLog.gammahalf(1) == 1.772453850905516
    @test PolyLog.gammahalf(3) == 0.88622692545275801
    @test PolyLog.gammahalf(343) == 9.4833675668247993e307
    @test PolyLog.gammahalf(345) == Inf
end
