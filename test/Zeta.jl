@testset "zeta" begin
    @test PolyLog.zeta(-263) == Inf
    @test PolyLog.zeta(-262) == 0.0
    @test PolyLog.zeta(-261) == -Inf
    @test PolyLog.zeta(-259) == 8.7601563446229215e306
    @test PolyLog.zeta(-257) == -5.1754977470366798e303
    @test PolyLog.zeta(-255) == 3.1055517596048927e300
    @test PolyLog.zeta(-4) == 0.0
    @test PolyLog.zeta(-3) == 1/120
    @test PolyLog.zeta(-2) == 0.0
    @test PolyLog.zeta(-1) == -1/12
    @test PolyLog.zeta(0) == -0.5
    @test PolyLog.zeta(1) == Inf
    @test PolyLog.zeta(2) == pi^2/6
    @test PolyLog.zeta(32) == 1.0000000002328312
    @test PolyLog.zeta(33) == 1.0000000001164155
    @test PolyLog.zeta(34) == 1.0000000000582077
    @test PolyLog.zeta(35) == 1.0000000000291039
end

@testset "zetahalf" begin
    @test PolyLog.zetahalf(-521) == Inf
    @test PolyLog.zetahalf(-519) == 3.982766118112322e307
    @test PolyLog.zetahalf(-3) == -0.025485201889833036
    @test PolyLog.zetahalf(-1) == -0.20788622497735457
    @test PolyLog.zetahalf(1) == -1.4603545088095868
    @test PolyLog.zetahalf(3) == 2.6123753486854883
    @test PolyLog.zetahalf(107) == 1.0000000000000001
    @test PolyLog.zetahalf(109) == 1.0
end
