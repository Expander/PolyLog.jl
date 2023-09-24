@testset "fac" begin
    @test_throws DomainError PolyLog.fac(-1)
    @test PolyLog.fac(  0) == 1.0
    @test PolyLog.fac(  1) == 1.0
    @test PolyLog.fac(  2) == 2.0
    @test PolyLog.fac(  3) == 6.0
    @test PolyLog.fac(  4) == 24.0
    @test PolyLog.fac(170) == 7.2574156153079990e306
    @test PolyLog.fac(171) == Inf

    @test PolyLog.fac(  0, BigFloat) == 1.0
    @test PolyLog.fac(171, BigFloat) == BigFloat(factorial(big(171)))
end

@testset "inv_fac" begin
    @test_throws DomainError PolyLog.inv_fac(-1)
    @test PolyLog.inv_fac(  0) == 1.0
    @test PolyLog.inv_fac(  1) == 1.0
    @test PolyLog.inv_fac(  2) == 1/2
    @test PolyLog.inv_fac(  3) == 1/6
    @test PolyLog.inv_fac(  4) == 1/24
    @test PolyLog.inv_fac(176) == 5.0529776809556651e-321
    @test PolyLog.inv_fac(177) == 2.8547896502574379e-323
    @test PolyLog.inv_fac(178) == 0.0

    @test PolyLog.inv_fac(178, BigFloat) == BigFloat(inv(factorial(big(178))))
end
