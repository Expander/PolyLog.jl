# zeta[n] for n = 2,...,21
const zetas = (
    zeta2, zeta3, zeta4, zeta5, zeta6,
    1.0083492773819228, 1.0040773561979443, 1.0020083928260822,
    1.0009945751278181, 1.0004941886041195, 1.0002460865533080,
    1.0001227133475785, 1.0000612481350587, 1.0000305882363070,
    1.0000152822594087, 1.0000076371976379, 1.0000038172932650,
    1.0000019082127166, 1.0000009539620339, 1.0000004769329868
)

function zeta(n::Integer)::Float64
    n <= 1 && throw(DomainError(n, "zeta(n) undefined for n <= 1"))
    n - 1 <= length(zetas) && return zetas[n - 1]

    sum = one(Float64)

    for k in 1:typemax(n)
        old_sum = sum
        sum += inv(2*k + 1)^n
        sum == old_sum && break
    end

    sum/(one(Float64) - 2.0^(-n))
end
