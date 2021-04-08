using Test
using DelimitedFiles

include("Common.jl")
include("../src/Li4.jl")

li4_data = open(readdlm, "data/Li4.txt")

for r in 1:size(li4_data, 1)
    row = li4_data[r, :]

    z            = row[1] + row[2]*1im
    li4_expected = row[3] + row[4]*1im

    @test is_equal(li4(z), li4_expected, 1e-14)
end
