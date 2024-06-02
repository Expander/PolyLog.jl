module PolyLogForwardDiffExt

using PolyLog
using ForwardDiff
using ForwardDiff: Dual, partials
function PolyLog.reli4(d::Dual{T}) where T
    val = ForwardDiff.value(d)
    x = reli4(val)
    dx = reli3(val)/val
    return Dual{T}(x, dx * partials(d))
end

function PolyLog.reli3(d::Dual{T}) where T
    val = ForwardDiff.value(d)
    x = reli3(val)
    dx = reli2(val)/val
    return Dual{T}(x, dx * partials(d))
end

function PolyLog.reli2(d::Dual{T}) where T
    val = ForwardDiff.value(d)
    x = reli2(val)
    dx = reli1(val)/val
    return Dual{T}(x, dx * partials(d))
end

function PolyLog.reli1(d::Dual{T}) where T
    val = ForwardDiff.value(d)
    x = reli1(val)
    dx = one(val) / (one(val) - val)
    return Dual{T}(x, dx*partials(d))
end

function PolyLog.reli(n::Integer,d::Dual{T}) where T
    val = ForwardDiff.value(d)
    x = reli(val)
    dx = PolyLog.reli(n-1,val)/val
    return Dual{T}(x, dx*partials(d))
end

end #module