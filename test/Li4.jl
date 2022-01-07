@testset "li4" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Li4.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        z        = row[1] + row[2]*1im
        expected = row[3] + row[4]*1im

        @test PolyLog.li4(z) ≈ expected atol=1e-14
    end
end
