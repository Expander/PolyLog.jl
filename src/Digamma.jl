# digamma for integer n > 0
function digamma(n::Integer)::Float64
    # Table[BernoulliB[2n]/(2 n), {n,1,8}]
    C = (
        0.083333333333333333, -0.0083333333333333333, 0.0039682539682539683, 
        -0.0041666666666666667, 0.0075757575757575758, -0.021092796092796093, 
        0.083333333333333333, -0.44325980392156863
    )
    n <= 0 && throw(DomainError(n, "digamma not implemented for n <= 0"))
    res = zero(Float64)
    if n < 7 # recurrence formula
        k = 7 - n
        for nu = 1:(k - 1)
            res -= inv(n + nu)
        end
        res -= inv(n)
        n += k
    end
    t = inv(n)
    res += log(n) - 0.5*t
    t *= t # 1/n^2
    res -= t*(C[1] + t*(C[2] + t*(C[3] + t*(C[4] + t*(C[5] + t*(C[6] + t*(C[7] + t*C[8])))))))
end
