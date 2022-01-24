PolyLog.jl
==========

[![test](https://github.com/Expander/PolyLog.jl/actions/workflows/build.yml/badge.svg)](https://github.com/Expander/PolyLog.jl/actions/workflows/build.yml)

The PolyLog.jl package provides Julia implementations of real and
complex polylogarithms, including the real and complex dilogarithm and
trilogarithm.


Example
-------

```.jl
using PolyLog

# real polylogarithms
li0(1.0)          # Re[Li_0(x)]
li1(1.0)          # Re[Li_1(x)]
li2(1.0)          # Re[Li_2(x)] (dilogarithm)
li3(1.0)          # Re[Li_3(x)] (trilogarithm)
li4(1.0)          # Re[Li_4(x)]
li(10, 1.0)       # Re[Li_n(x)] for all integers n >= 0 (here: n = 10)

# complex polylogarithms
li0(1.0 + 1.0im)  # Li_0(z)
li1(1.0 + 1.0im)  # Li_1(z)
li2(1.0 + 1.0im)  # Li_2(z) (dilogarithm)
li3(1.0 + 1.0im)  # Li_3(z) (trilogarithm)
li4(1.0 + 1.0im)  # Li_4(z)
li5(1.0 + 1.0im)  # Li_5(z)
li6(1.0 + 1.0im)  # Li_6(z)
```


Documentation
-------------

[https://docs.juliahub.com/PolyLog/](https://docs.juliahub.com/PolyLog/)


Notes
-----

The implementation of the complex dilogarithm has been inspired by the
implementation in [SPheno](https://spheno.hepforge.org) and has been
translated to Julia.


Citation
--------

~~~.bibtex
@software{PolyLog.jl,
    author       = {{Alexander Voigt}},
    title        = {{PolyLog.jl}},
    year         = {2021},
    version      = {1.7.0},
    url          = {https://github.com/Expander/PolyLog.jl},
    note         = {[License: MIT]}
}
~~~


Copying
-------

PolyLog.jl is licenced under the MIT License.


Links
-----

Refer to the package
[Polylogarithms.jl](https://github.com/mroughan/Polylogarithms.jl) for
a Julia implementation of polylogarithms of arbitrary order.
