@testset "ln_sqr" begin
    for x in -10:10
        @test PolyLog.ln_sqr(Float64(x)) == abs2(log(Float64(x) + 0.0im))
    end
end
