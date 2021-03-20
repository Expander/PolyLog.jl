include("../src/Li2.jl")

n = 10_000_000
z_min = -5
z_max = 5
cmpl_data = (z_max - z_min)*rand(ComplexF64, n) + z_min*(1.0 + 1.0im)*ones(n)
real_data = map(real, cmpl_data)

println("Benchmarking li2::Float64")

map(li2, real_data)       # trigger compilation
@time map(li2, real_data)
@time map(li2, real_data)
@time map(li2, real_data)

println("Benchmarking li2::ComplexF64")

map(li2, cmpl_data)       # trigger compilation
@time map(li2, cmpl_data)
@time map(li2, cmpl_data)
@time map(li2, cmpl_data)
