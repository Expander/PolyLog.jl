@testset "li2" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Li2.txt"))

    for r in 1:size(data, 1)
        row = data[r, :]

        z            = row[1] + row[2]*1im
        li2_expected = row[3] + row[4]*1im

        if imag(z) == 0.0
            @test is_equal(PolyLog.li2(real(z)), real(li2_expected), 1e-14)
        end

        @test is_equal(PolyLog.li2(z), li2_expected, 1e-14)
    end

    @test PolyLog.li2(1.0) == 1.6449340668482264
end
