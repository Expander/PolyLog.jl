import PolyLog

n = 10_000_000
z_min = -5
z_max = 5
c64_data = (z_max - z_min)*rand(ComplexF64, n) + z_min*(1.0 + 1.0im)*ones(n)
f64_data = map(real, c64_data)
f32_data = map(Float32, f64_data)

time_reli2(data) = @time map(PolyLog.reli2, data)
time_reli3(data) = @time map(PolyLog.reli3, data)
time_reli4(data) = @time map(PolyLog.reli4, data)
time_reli5(data) = @time map(PolyLog.reli5, data)
time_reli6(data) = @time map(PolyLog.reli6, data)

time_li2(data) = @time map(PolyLog.li2, data)
time_li3(data) = @time map(PolyLog.li3, data)
time_li4(data) = @time map(PolyLog.li4, data)
time_li5(data) = @time map(PolyLog.li5, data)
time_li6(data) = @time map(PolyLog.li6, data)

println("Benchmarking li2::Float64")

map(PolyLog.reli2, f64_data)       # trigger compilation
time_reli2(f64_data)
time_reli2(f64_data)
time_reli2(f64_data)

println("Benchmarking li2::Float32")

map(PolyLog.reli2, f32_data)       # trigger compilation
time_reli2(f32_data)
time_reli2(f32_data)
time_reli2(f32_data)

println("Benchmarking li2::ComplexF64")

map(PolyLog.li2, c64_data)       # trigger compilation
time_li2(c64_data)
time_li2(c64_data)
time_li2(c64_data)

println("Benchmarking li3::Float64")

map(PolyLog.reli3, f64_data)       # trigger compilation
time_reli3(f64_data)
time_reli3(f64_data)
time_reli3(f64_data)

println("Benchmarking li3::ComplexF64")

map(PolyLog.li3, c64_data)       # trigger compilation
time_li3(c64_data)
time_li3(c64_data)
time_li3(c64_data)

println("Benchmarking reli4::Float64")

map(PolyLog.reli4, f64_data)       # trigger compilation
time_reli4(f64_data)
time_reli4(f64_data)
time_reli4(f64_data)

println("Benchmarking li4::ComplexF64")

map(PolyLog.li4, c64_data)       # trigger compilation
time_li4(c64_data)
time_li4(c64_data)
time_li4(c64_data)

println("Benchmarking li5::ComplexF64")

map(PolyLog.li5, c64_data)       # trigger compilation
time_li5(c64_data)
time_li5(c64_data)
time_li5(c64_data)

println("Benchmarking li6::ComplexF64")

map(PolyLog.li6, c64_data)       # trigger compilation
time_li6(c64_data)
time_li6(c64_data)
time_li6(c64_data)

println("Benchmarking li::Float64")

time_reli(k, data) = @time map(x -> PolyLog.reli(k, x), data)

n = 1_000_000
c64_data = (z_max - z_min)*rand(ComplexF64, n) + z_min*(1.0 + 1.0im)*ones(n)
f64_data = map(real, c64_data)

for k in vcat(collect(-10:2:10), [100, 1000, 1000_000, -100, -1000, -1000_000])
    println("Benchmarking li($(k),x)::Float64")
    map(x -> PolyLog.reli(k, x), f64_data) # trigger compilation
    time_reli(k, f64_data)                 # trigger compilation
    time_reli(k, f64_data)
    time_reli(k, f64_data)
end

println("Benchmarking li::ComplexF64")

time_li(k, data) = @time map(x -> PolyLog.li(k, x), data)

for k in vcat(collect(-10:2:10), [100, 1000, 1000_000, -100, -1000, -1000_000])
    println("Benchmarking li($(k),x)::ComplexF64")
    map(x -> PolyLog.li(k, x), c64_data) # trigger compilation
    time_li(k, c64_data)                 # trigger compilation
    time_li(k, c64_data)
    time_li(k, c64_data)
end
