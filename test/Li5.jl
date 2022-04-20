@testset "li5" begin
    data = read_from(joinpath(@__DIR__, "data", "Li5.txt"))

    for i in 1:size(data, 1)
        z = data[i,1]
        expected = data[i,2]

        @test PolyLog.li5(z) ≈ expected atol=1e-14
    end
end
