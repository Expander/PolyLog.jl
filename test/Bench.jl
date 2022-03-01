import PolyLog

n = 10_000_000
z_min = -5
z_max = 5
cmpl_data = (z_max - z_min)*rand(ComplexF64, n) + z_min*(1.0 + 1.0im)*ones(n)
real_data = map(real, cmpl_data)

time_li2(data) = @time map(PolyLog.li2, data)
time_li3(data) = @time map(PolyLog.li3, data)
time_li4(data) = @time map(PolyLog.li4, data)
time_li5(data) = @time map(PolyLog.li5, data)
time_li6(data) = @time map(PolyLog.li6, data)

println("Benchmarking li2::Float64")

map(PolyLog.li2, real_data)       # trigger compilation
time_li2(real_data)
time_li2(real_data)
time_li2(real_data)

println("Benchmarking li2::ComplexF64")

map(PolyLog.li2, cmpl_data)       # trigger compilation
time_li2(cmpl_data)
time_li2(cmpl_data)
time_li2(cmpl_data)

println("Benchmarking li3::Float64")

map(PolyLog.li3, real_data)       # trigger compilation
time_li3(real_data)
time_li3(real_data)
time_li3(real_data)

println("Benchmarking li3::ComplexF64")

map(PolyLog.li3, cmpl_data)       # trigger compilation
time_li3(cmpl_data)
time_li3(cmpl_data)
time_li3(cmpl_data)

println("Benchmarking li4::Float64")

map(PolyLog.li4, real_data)       # trigger compilation
time_li4(real_data)
time_li4(real_data)
time_li4(real_data)

println("Benchmarking li4::ComplexF64")

map(PolyLog.li4, cmpl_data)       # trigger compilation
time_li4(cmpl_data)
time_li4(cmpl_data)
time_li4(cmpl_data)

println("Benchmarking li5::ComplexF64")

map(PolyLog.li5, cmpl_data)       # trigger compilation
time_li5(cmpl_data)
time_li5(cmpl_data)
time_li5(cmpl_data)

println("Benchmarking li6::ComplexF64")

map(PolyLog.li6, cmpl_data)       # trigger compilation
time_li6(cmpl_data)
time_li6(cmpl_data)
time_li6(cmpl_data)

println("Benchmarking li::Float64")

time_li(k, data) = @time map(x -> PolyLog.li(k, x), data)

n = 1_000_000
cmpl_data = (z_max - z_min)*rand(ComplexF64, n) + z_min*(1.0 + 1.0im)*ones(n)
real_data = map(real, cmpl_data)

for k in vcat(collect(-10:2:10), [100, 1000, 1000_000, -100, -1000, -1000_000])
    println("Benchmarking li($(k),x)::Float64")
    map(x -> PolyLog.li(k, x), real_data) # trigger compilation
    time_li(k, real_data)                 # trigger compilation
    time_li(k, real_data)
    time_li(k, real_data)
end

println("Benchmarking li::ComplexF64")

for k in vcat(collect(-10:2:10), [100, 1000, 1000_000, -100, -1000, -1000_000])
    println("Benchmarking li($(k),x)::ComplexF64")
    map(x -> PolyLog.li(k, x), cmpl_data) # trigger compilation
    time_li(k, cmpl_data)                 # trigger compilation
    time_li(k, cmpl_data)
    time_li(k, cmpl_data)
end
