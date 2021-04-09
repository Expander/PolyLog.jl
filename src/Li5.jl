"""
    li5(z::ComplexF64)::ComplexF64

Returns the complex 5th order polylogarithm of a complex number `z` of type `ComplexF64`.

Author: Alexander Voigt

License: MIT

# Example
```julia
li5(1.0 + 1.0im)
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
            return z4 + 0.0im
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
        v4::ComplexF64 = v2*v2
        v8::ComplexF64 = v4*v4
        c1::Float64 = 1.2020569031595943 # zeta(3)
        c2::Float64 = 0.82246703342411322
        c3::ComplexF64 = (11.0/6.0 - clog(-v))/6.0
        c4::Float64 = -1.0/48.0

        cs = (
            -6.9444444444444444e-04, 1.6534391534391534e-06,
            -1.0935444136502338e-08, 1.0438378493934049e-10,
            -1.2165942300622435e-12, 1.6130006528350101e-14,
            -2.3428810452879340e-16
        )

        return z4 + v2 * (c2 + v2 * c4) +
            v * (
                c1 +
                c3*v2 +
                v4*(cs[1] + v2*cs[2]) +
                v8*(cs[3] + v2*cs[4] + v4*(cs[5] + v2*cs[6])) +
                v8*v8*cs[7]
            )
    end

    (u::ComplexF64, rest::ComplexF64, sgn::Float64) = if nz <= 1.0
        (-clog(1.0 - z), 0.0 + 0.0im, 1.0)
    else # nz > 1.0
        arg::Float64 = pz > 0.0 ? pz - pi : pz + pi
        lmz::ComplexF64 = lnz + arg*im # clog(z)
        lmz2::ComplexF64 = lmz*lmz
        (-clog(1.0 - 1.0/z), 1.0/360.0*(-7.0*pi^4 + lmz2*(-30.0*pi^2 - 15.0*lmz2)), -1.0)
    end

    u2::ComplexF64 = u*u
    u4::ComplexF64 = u2*u2
    u8::ComplexF64 = u4*u4

    rest + sgn * (
       u*bf[1] +
       u2*(bf[2] + u*bf[3]) +
       u4*(bf[4] + u*bf[5] + u2*(bf[6] + u*bf[7])) +
       u8*(bf[8] + u*bf[9] + u2*(bf[10] + u*bf[11]) +
           u4*(bf[12] + u*bf[13] + u2*(bf[14] + u*bf[15]))) +
       u8*u8*(bf[16] + u*bf[17] + u2*bf[18])
    )
end
