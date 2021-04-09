li3_data = open(readdlm, "data/Li3.txt")

for r in 1:size(li3_data, 1)
    row = li3_data[r, :]

    z            = row[1] + row[2]*1im
    li3_expected = row[3] + row[4]*1im

    @test is_equal(PolyLog.li3(z), li3_expected, 1e-14)
end
