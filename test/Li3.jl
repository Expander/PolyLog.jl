@testset "li3" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Li3.txt"))

    for r in 1:size(data, 1)
        row = data[r, :]

        z            = row[1] + row[2]*1im
        li3_expected = row[3] + row[4]*1im

        @test is_equal(PolyLog.li3(z), li3_expected, 1e-14)
    end
end
