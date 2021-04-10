using PolyLog

n = 10_000_000
z_min = -5
z_max = 5
cmpl_data = (z_max - z_min)*rand(ComplexF64, n) + z_min*(1.0 + 1.0im)*ones(n)
real_data = map(real, cmpl_data)

println("Benchmarking li2::Float64")

map(PolyLog.li2, real_data)       # trigger compilation
@time map(PolyLog.li2, real_data)
@time map(PolyLog.li2, real_data)
@time map(PolyLog.li2, real_data)

println("Benchmarking li2::ComplexF64")

map(PolyLog.li2, cmpl_data)       # trigger compilation
@time map(PolyLog.li2, cmpl_data)
@time map(PolyLog.li2, cmpl_data)
@time map(PolyLog.li2, cmpl_data)

println("Benchmarking li3::ComplexF64")

map(PolyLog.li3, cmpl_data)       # trigger compilation
@time map(PolyLog.li3, cmpl_data)
@time map(PolyLog.li3, cmpl_data)
@time map(PolyLog.li3, cmpl_data)

println("Benchmarking li4::ComplexF64")

map(PolyLog.li4, cmpl_data)       # trigger compilation
@time map(PolyLog.li4, cmpl_data)
@time map(PolyLog.li4, cmpl_data)
@time map(PolyLog.li4, cmpl_data)

println("Benchmarking li5::ComplexF64")

map(PolyLog.li5, cmpl_data)       # trigger compilation
@time map(PolyLog.li5, cmpl_data)
@time map(PolyLog.li5, cmpl_data)
@time map(PolyLog.li5, cmpl_data)

println("Benchmarking li6::ComplexF64")

map(PolyLog.li6, cmpl_data)       # trigger compilation
@time map(PolyLog.li6, cmpl_data)
@time map(PolyLog.li6, cmpl_data)
@time map(PolyLog.li6, cmpl_data)
