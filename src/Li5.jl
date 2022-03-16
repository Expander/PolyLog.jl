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
    if imag(z) == 0.0
        if real(z) == 0.0
            return 0.0 + 0.0im
        end
        if real(z) == 1.0
            return zeta5 + 0.0im
        end
        if real(z) == -1.0
            return -15/16*zeta5 + 0.0im
        end
    end

    nz::Float64  = abs2(z)
    pz::Float64  = angle(z)
    lnz::Float64 = 0.5*log(nz)

    if lnz*lnz + pz*pz < 1.0 # |log(z)| < 1
        cs = (
            -1.1574074074074074e-04, 2.0667989417989418e-07,
            -1.0935444136502338e-09, 8.6986487449450412e-12,
            -8.6899587861588824e-14, 1.0081254080218813e-15
        )

        u  = lnz + pz*im
        u2 = u*u
        c1 = zeta4
        c2 = 0.5*zeta3
        c3 = 0.27415567780803774
        c4 = 1/24*(25/12 - clog(-u))
        c5 = -1/240

        zeta5 + u*c1 +
        u2*(c2 + u*c3 +
        u2*(c4 + u*c5 +
        u2*(cs[1] +
        u2*(cs[2] +
        u2*(cs[3] +
        u2*(cs[4] +
        u2*(cs[5] +
        u2*(cs[6]))))))))
    else
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

        (u, rest) = if nz <= 1.0
            (-clog(1.0 - z), 0.0 + 0.0im)
        else # nz > 1.0
            arg = pz > 0.0 ? pz - pi : pz + pi
            lmz = lnz + arg*im # clog(z)
            lmz2 = lmz*lmz
            (-clog(1.0 - inv(z)), lmz*(-7/4*zeta4 - lmz2*(1/6*zeta2 + 1/120*lmz2)))
        end

        u2 = u*u
        u4 = u2*u2
        u8 = u4*u4

        rest +
        u*bf[1] +
        u2*(bf[2] + u*bf[3]) +
        u4*(bf[4] + u*bf[5] + u2*(bf[6] + u*bf[7])) +
        u8*(bf[8] + u*bf[9] + u2*(bf[10] + u*bf[11]) +
            u4*(bf[12] + u*bf[13] + u2*(bf[14] + u*bf[15]))) +
        u8*u8*(bf[16] + u*bf[17] + u2*(bf[18] + u*bf[19]))
    end
end
