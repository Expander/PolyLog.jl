using Test
using DelimitedFiles

include("Common.jl")
include("../src/Li2.jl")

li2_data = open(readdlm, "test/data/Li2.txt")

for r in 1:size(li2_data, 1)
    row = li2_data[r, :]

    z            = row[1] + row[2]*1im
    li2_expected = row[3] + row[4]*1im

    if imag(z) == 0.0
        # println("z = ", z, " li2_expected = ", real(li2_expected), " li2_have = ", li2(real(z)))
        @test is_equal(li2(real(z)), real(li2_expected), 1e-14)
    end
end

@test li2(1.0) == 1.6449340668482264
