@testset "li5" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Li5.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        z        = row[1] + row[2]*1im
        expected = row[3] + row[4]*1im

        @test PolyLog.li5(z) ≈ expected atol=1e-14
    end
end
