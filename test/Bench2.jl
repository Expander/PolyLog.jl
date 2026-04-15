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

benchmark_and_print(PolyLog.li1, 1.1        , "li1::Float64")
benchmark_and_print(PolyLog.li1, 1.1 + 0.5im, "li1::ComplexF64")
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

const MAX_DECIMAL_DIGITS = 40
const MAX_BINARY_DIGITS = ceil(Integer, MAX_DECIMAL_DIGITS*log(10)/log(2))

setprecision(BigFloat, ceil(Integer, MAX_BINARY_DIGITS)) do
    for n in [-1000_000, -1000, -100, -10, -2, -1, 0, 1, 2, 3, 4, 5, 6, 10, 100, 1000, 1000_000]
        for z in [BigFloat("0.5") + BigFloat("0.5")*1im,
                  BigFloat("1.1") + BigFloat("1.1")*1im]
            x = real(z)
            benchmark_and_print(x -> PolyLog.li(n, x), x, "li($(n), $(x))::BigFloat")
            benchmark_and_print(z -> PolyLog.li(n, z), z, "li($(n), $(z))::Complex{BigFloat}")
        end
    end
end
