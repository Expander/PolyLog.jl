@testset "Missing" begin
    @test ismissing(PolyLog.li0(missing))
    @test ismissing(PolyLog.li1(missing))
    @test ismissing(PolyLog.li2(missing))
    @test ismissing(PolyLog.li3(missing))
    @test ismissing(PolyLog.li4(missing))
    @test ismissing(PolyLog.li5(missing))
    @test ismissing(PolyLog.li6(missing))

    @test ismissing(PolyLog.reli1(missing))
    @test ismissing(PolyLog.reli2(missing))
    @test ismissing(PolyLog.reli3(missing))
    @test ismissing(PolyLog.reli4(missing))
    
    for n in -16:16
        @test ismissing(PolyLog.li(n, missing))
        @test ismissing(PolyLog.reli(n, missing))
    end
end
