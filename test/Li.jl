@testset "li" begin
    for n in vcat(collect(0:6), [100])
        data = open(readdlm, joinpath(@__DIR__, "data", "Li$(n).txt"))

        for r in 1:size(data, 1)
            row      = data[r, :]
            z        = row[1] + row[2]*1im
            expected = row[3] + row[4]*1im

            if imag(z) == 0.0
                @test PolyLog.li(n, real(z)) ≈ real(expected) atol=1e-14
            end
        end
    end

    @test_throws DomainError PolyLog.li(-1, 1.0)
end
