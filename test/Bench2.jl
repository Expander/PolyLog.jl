using BenchmarkTools
import PolyLog

function benchmark_and_print(f, z, str)
    println("------------------------------------------")
    println("Benchmarking ", str)
    println("------------------------------------------")
    b = @benchmark $f($z)
    io = IOBuffer()
    show(io, "text/plain", b)
    println(String(take!(io)))
end

benchmark_and_print(PolyLog.li2, 1.1        , "li2::Float64")
benchmark_and_print(PolyLog.li2, 1.1 + 0.5im, "li2::ComplexF64")
benchmark_and_print(PolyLog.li3, 1.1        , "li3::Float64")
benchmark_and_print(PolyLog.li3, 1.1 + 0.5im, "li3::ComplexF64")
benchmark_and_print(PolyLog.li4, 1.1        , "li4::Float64")
benchmark_and_print(PolyLog.li4, 1.1 + 0.5im, "li4::ComplexF64")
benchmark_and_print(PolyLog.li5, 1.1 + 0.5im, "li5::ComplexF64")
benchmark_and_print(PolyLog.li6, 1.1 + 0.5im, "li6::ComplexF64")

for n in vcat(collect(2:4), [100, 1000, 1000_000])
    benchmark_and_print(x -> PolyLog.li(n, x), -1.5, "li($(n), -1.5)::Float64")
    benchmark_and_print(x -> PolyLog.li(n, x),  0.4, "li($(n), 0.4)::Float64")
    benchmark_and_print(x -> PolyLog.li(n, x),  1.5, "li($(n), 1.5)::Float64")
end
