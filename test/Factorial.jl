@testset "fac" begin
    @test_throws DomainError PolyLog.fac(-1, Float64)
    @test_throws DomainError PolyLog.fac(-1, BigFloat)

    @test typeof(PolyLog.fac(1, Float16)) == Float16
    @test typeof(PolyLog.fac(1, Float32)) == Float32
    @test typeof(PolyLog.fac(1, Float64)) == Float64
    @test typeof(PolyLog.fac(1, BigFloat)) == BigFloat

    @test PolyLog.fac(  0, Float64) == 1.0
    @test PolyLog.fac(  1, Float64) == 1.0
    @test PolyLog.fac(  2, Float64) == 2.0
    @test PolyLog.fac(  3, Float64) == 6.0
    @test PolyLog.fac(  4, Float64) == 24.0
    @test PolyLog.fac(170, Float64) == 7.2574156153079990e306
    @test PolyLog.fac(171, Float64) == Inf

    @test PolyLog.fac(  0, BigFloat) == BigFloat("1.0")
    @test PolyLog.fac(  1, BigFloat) == BigFloat("1.0")
    @test PolyLog.fac(  2, BigFloat) == BigFloat("2.0")
    @test PolyLog.fac(  3, BigFloat) == BigFloat("6.0")
    @test PolyLog.fac(  4, BigFloat) == BigFloat("24.0")
    @test PolyLog.fac(171, BigFloat) == BigFloat(factorial(big(171)))
end

@testset "inv_fac" begin
    @test_throws DomainError PolyLog.inv_fac(-1, Float64)
    @test_throws DomainError PolyLog.inv_fac(-1, BigFloat)

    @test typeof(PolyLog.inv_fac(1, Float16)) == Float16
    @test typeof(PolyLog.inv_fac(1, Float32)) == Float32
    @test typeof(PolyLog.inv_fac(1, Float64)) == Float64
    @test typeof(PolyLog.inv_fac(1, BigFloat)) == BigFloat

    @test PolyLog.inv_fac(  0, Float64) == 1.0
    @test PolyLog.inv_fac(  1, Float64) == 1.0
    @test PolyLog.inv_fac(  2, Float64) == 1/2
    @test PolyLog.inv_fac(  3, Float64) == 1/6
    @test PolyLog.inv_fac(  4, Float64) == 1/24
    @test PolyLog.inv_fac(176, Float64) == 5.0529776809556651e-321
    @test PolyLog.inv_fac(177, Float64) == 2.8547896502574379e-323
    @test PolyLog.inv_fac(178, Float64) == 0.0

    @test PolyLog.inv_fac(  0, BigFloat) == BigFloat("1.0")
    @test PolyLog.inv_fac(  1, BigFloat) == BigFloat("1.0")
    @test PolyLog.inv_fac(  2, BigFloat) == BigFloat("1.0")/2
    @test PolyLog.inv_fac(  3, BigFloat) == BigFloat("1.0")/6
    @test PolyLog.inv_fac(  4, BigFloat) == BigFloat("1.0")/24
    @test PolyLog.inv_fac(178, BigFloat) == BigFloat(inv(factorial(big(178))))
end
