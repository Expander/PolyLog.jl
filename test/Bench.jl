include("../src/Li2.jl")

n = 1000_000
z_min = -5
z_max = 5
data = (z_max - z_min)*rand(ComplexF64, n) + z_min*(1.0 + 1.0im)*ones(n)

@time map(li2, data)
