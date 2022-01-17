"""
    li5(z::ComplexF64)::ComplexF64

Returns the complex 5th order polylogarithm
``\\operatorname{Li}_5(z)`` of a complex number ``z`` of type
`ComplexF64`.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using PolyLog)
julia> li5(1.0 + 1.0im)
0.9874666591701124 + 1.068441607107422im
```
"""
function li5(z::ComplexF64)::ComplexF64
    function clog(z)
        az::Float64 = angle(z)
        return 0.5*log(abs2(z)) + (imag(z) == 0.0 && az < 0.0 ? -az : az)*1.0im
    end

    z5::Float64 = 1.0369277551433699
    bf = (
        1.0                   , -15.0/32.0             ,
        1.3953189300411523e-01, -2.8633777006172840e-02,
        4.0317412551440329e-03, -3.3985018004115226e-04,
        4.5445184621617666e-06,  2.3916808048569012e-06,
       -1.2762692600122747e-07, -3.1628984306505932e-08,
        3.2848118445335192e-09,  4.7613713995660579e-10,
       -8.0846898171909830e-11, -7.2387648587737207e-12,
        1.9439760115173968e-12,  1.0256978405977236e-13,
       -4.6180551009884830e-14, -1.1535857196470580e-15,
        1.0903545401333394e-15
    )

    if imag(z) == 0.0
        if real(z) == 0.0
            return 0.0 + 0.0im
        end
        if real(z) == 1.0
            return z5 + 0.0im
        end
        if real(z) == -1.0
            return -15.0/16.0*z5 + 0.0im
        end
    end

    nz::Float64  = abs2(z)
    pz::Float64  = angle(z)
    lnz::Float64 = 0.5*log(nz)

    if lnz*lnz + pz*pz < 1.0 # |log(z)| < 1
        v::ComplexF64  = lnz + pz*im
        v2::ComplexF64 = v*v
        c1::Float64 = 1.0823232337111382 # zeta(4)
        c2::Float64 = 0.60102845157979714 # zeta(3)/2
        c3::Float64 = 0.27415567780803774
        c4::ComplexF64 = (25.0/12.0 - clog(-v))/24.0
        c5::Float64 = -1.0/240.0

        cs = (
            -1.1574074074074074e-04, 2.0667989417989418e-07,
            -1.0935444136502338e-09, 8.6986487449450412e-12,
            -8.6899587861588824e-14, 1.0081254080218813e-15
        )

        return z5 + v * c1 +
            v2 * (c2 + v * c3 +
            v2 * (c4 + v * c5 +
            v2 * (cs[1] +
            v2 * (cs[2] +
            v2 * (cs[3] +
            v2 * (cs[4] +
            v2 * (cs[5] +
            v2 * (cs[6]))))))))
    end

    (u::ComplexF64, rest::ComplexF64) = if nz <= 1.0
        (-clog(1.0 - z), 0.0 + 0.0im)
    else # nz > 1.0
        arg::Float64 = pz > 0.0 ? pz - pi : pz + pi
        lmz::ComplexF64 = lnz + arg*im # clog(z)
        lmz2::ComplexF64 = lmz*lmz
        (-clog(1.0 - 1.0/z), -1.0/360.0*lmz*(7.0*pi^4 + lmz2*(10.0*pi^2 + 3.0*lmz2)))
    end

    u2::ComplexF64 = u*u
    u4::ComplexF64 = u2*u2
    u8::ComplexF64 = u4*u4

    rest +
    u*bf[1] +
    u2*(bf[2] + u*bf[3]) +
    u4*(bf[4] + u*bf[5] + u2*(bf[6] + u*bf[7])) +
    u8*(bf[8] + u*bf[9] + u2*(bf[10] + u*bf[11]) +
        u4*(bf[12] + u*bf[13] + u2*(bf[14] + u*bf[15]))) +
    u8*u8*(bf[16] + u*bf[17] + u2*(bf[18] + u*bf[19]))
end
