if isdefined(Base, :get_extension)
    import ForwardDiff
    import ChainRulesTestUtils

    @testset "ForwardDiff" begin
        @test ForwardDiff.derivative(PolyLog.li0, float(pi)) == 1/(1 - pi)^2
        @test ForwardDiff.derivative(PolyLog.li0, 0.0) == 1
        @test ForwardDiff.derivative(PolyLog.li0, 1.0) == Inf

        @test ForwardDiff.derivative(PolyLog.reli1, float(pi)) == 1/(1 - pi)
        @test ForwardDiff.derivative(PolyLog.reli1, 0.0) == 1.0
        ChainRulesTestUtils.test_frule(PolyLog.reli1, 0.0)
        ChainRulesTestUtils.test_rrule(PolyLog.reli1, float(pi))

        @test ForwardDiff.derivative(PolyLog.reli2, float(pi)) == PolyLog.reli1(pi)/pi
        @test ForwardDiff.derivative(PolyLog.reli2, 0.0) == 1.0
        ChainRulesTestUtils.test_frule(PolyLog.reli2, 0.0)
        ChainRulesTestUtils.test_rrule(PolyLog.reli2, float(pi))

        @test ForwardDiff.derivative(PolyLog.reli3, float(pi)) == PolyLog.reli2(pi)/pi
        @test ForwardDiff.derivative(PolyLog.reli3, 0.0) == 1.0
        ChainRulesTestUtils.test_frule(PolyLog.reli3, 0.0)
        ChainRulesTestUtils.test_rrule(PolyLog.reli3, float(pi))

        @test ForwardDiff.derivative(PolyLog.reli4, float(pi)) == PolyLog.reli3(pi)/pi
        @test ForwardDiff.derivative(PolyLog.reli4, 0.0) == 1.0
        ChainRulesTestUtils.test_frule(PolyLog.reli4, 0.0)
        ChainRulesTestUtils.test_rrule(PolyLog.reli4, float(pi))

        for n in vcat(collect(-10:10), [100, 1000000])
            @test ForwardDiff.derivative(z -> PolyLog.reli(n, z), float(pi)) == PolyLog.reli(n - 1, pi)/pi
            @test ForwardDiff.derivative(z -> PolyLog.reli(n, z), 0.0) == 1.0
            ChainRulesTestUtils.test_frule(PolyLog.reli, n, 0.0)
            ChainRulesTestUtils.test_rrule(PolyLog.reli, n, float(pi))
        end
    end
end
