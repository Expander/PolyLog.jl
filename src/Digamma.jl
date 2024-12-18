# digamma for integer n > 0, following
# [K.S. Kölbig: Programs for computing the logarithm of the gamma
# function, and the digamma function, for complex argument, Computer
# Physics Communications, Volume 4, Issue 2, 1972, Pages 221-226, ISSN
# 0010-4655, https://doi.org/10.1016/0010-4655(72)90012-4]
function digamma(n::Integer)::Float64
    # Table[BernoulliB[2n]/(2 n), {n,1,8}]
    C = (
        0.083333333333333333, -0.0083333333333333333,  0.0039682539682539683,
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
    t2 = t*t # 1/n^2
    res + log(n) - one(Float64)/2*t - t2*(C[1] + t2*(C[2] + t2*(C[3] + t2*(C[4] + t2*(C[5] + t2*(C[6] + t2*(C[7] + t2*C[8])))))))
end
