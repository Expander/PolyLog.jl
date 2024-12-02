using Test
import PolyLog
if isdefined(Base, :get_extension)
    import ForwardDiff
    import ChainRulesTestUtils
end

include("TestPrecision.jl")
include("DataReader.jl")
include("DataTester.jl")
include("Digamma.jl")
include("Eta.jl")
include("Factorial.jl")
include("Harmonic.jl")
include("Li0.jl")
include("Li1.jl")
include("Li2.jl")
include("Li3.jl")
include("Li4.jl")
include("Li5.jl")
include("Li6.jl")
include("Li.jl")
include("Log.jl")
include("Missing.jl")
include("Zeta.jl")
