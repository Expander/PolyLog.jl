module PolyLogChainRulesExt

using PolyLog
using ChainRulesCore

ChainRulesCore.@scalar_rule PolyLog.reli4(x) ifelse(iszero(x), one(x), PolyLog.reli3(x)/x)
ChainRulesCore.@scalar_rule PolyLog.reli3(x) ifelse(iszero(x), one(x), PolyLog.reli2(x)/x)
ChainRulesCore.@scalar_rule PolyLog.reli2(x) ifelse(iszero(x), one(x), PolyLog.reli1(x)/x)
ChainRulesCore.@scalar_rule PolyLog.reli1(x) one(x)/(one(x) - x)
ChainRulesCore.@scalar_rule(PolyLog.reli(n, x), (ChainRulesCore.NoTangent(), ifelse(iszero(x), one(x), PolyLog.reli(n - 1, x)/x)))

end #module
