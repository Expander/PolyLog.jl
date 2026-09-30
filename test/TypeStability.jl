@testset "TypeStability ($T)" for T in (Float16, Float32, Float64, BigFloat)
    @test @inferred(PolyLog.li0(T(1))) isa T
    @test @inferred(PolyLog.li0(Complex{T}(1))) isa Complex{T}

    @test @inferred(PolyLog.li1(T(1))) isa Complex{T}
    @test @inferred(PolyLog.li1(Complex{T}(1))) isa Complex{T}
    @test @inferred(PolyLog.reli1(T(1))) isa T

    @test @inferred(PolyLog.li2(T(1))) isa Complex{T}
    @test @inferred(PolyLog.li2(Complex{T}(1))) isa Complex{T}
    @test @inferred(PolyLog.reli2(T(1))) isa T

    if T != BigFloat
        @test @inferred(PolyLog.li3(T(1))) isa Complex{T}
        @test @inferred(PolyLog.li3(Complex{T}(1))) isa Complex{T}
        @test @inferred(PolyLog.reli3(T(1))) isa T

        @test @inferred(PolyLog.li4(T(1))) isa Complex{T}
        @test @inferred(PolyLog.li4(Complex{T}(1))) isa Complex{T}
        @test @inferred(PolyLog.reli4(T(1))) isa T

        @test @inferred(PolyLog.li5(T(1))) isa Complex{T}
        @test @inferred(PolyLog.li5(Complex{T}(1))) isa Complex{T}

        @test @inferred(PolyLog.li6(T(1))) isa Complex{T}
        @test @inferred(PolyLog.li6(Complex{T}(1))) isa Complex{T}
    end

    for n in -10:10
        @test @inferred(PolyLog.li(n,T(1))) isa Complex{T}
        @test @inferred(PolyLog.li(n,Complex{T}(1))) isa Complex{T}
        @test @inferred(PolyLog.reli(n,T(1))) isa T
    end
end
