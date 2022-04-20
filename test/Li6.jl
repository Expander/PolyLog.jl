@testset "li6" begin
    data = read_from(joinpath(@__DIR__, "data", "Li6.txt"))

    for d in eachrow(data)
        z = d[1]
        expected = d[2]

        @test PolyLog.li6(z) ≈ expected atol=1e-14
    end
end
