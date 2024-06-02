module PolyLogChainRulesExt

using PolyLog
using ChainRulesCore
@scalar_rule reli4(x) reli3(x)/x
@scalar_rule reli3(x) reli2(x)/x
@scalar_rule reli2(x) reli1(x)/x
@scalar_rule reli1(x) one(x)/(one(x)-x)
@scalar_rule reli(n,x) reli(n-1,x)/x
end #module