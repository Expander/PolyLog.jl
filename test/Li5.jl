@testset "li5" begin
    data = read_from(joinpath(@__DIR__, "data", "Li5.txt"))

    for d in eachrow(data)
        z = d[1]
        expected = d[2]

        @test PolyLog.li5(z) ≈ expected atol=1e-14
    end
end
