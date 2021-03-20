include("../src/Li2.jl")

n = 10_000_000
z_min = -5
z_max = 5
data = (z_max - z_min)*rand(ComplexF64, n) + z_min*(1.0 + 1.0im)*ones(n)

println("Benchmarking li2::ComplexF64")

li2(1.0 + 1.0im) # trigger compilation
@time map(li2, data)
@time map(li2, data)
@time map(li2, data)
