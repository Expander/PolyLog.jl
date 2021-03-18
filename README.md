PolyLog.jl
==========

The PolyLog.jl package provides Rust implementations of real and
complex polylogarithms.


Example
-------

```
include("src/PolyLog.jl")

println(Li2(2.0))
```


Notes
-----

The implementation of the complex dilogarithm has been inspired by the
implementation in [SPheno](https://spheno.hepforge.org) and has been
translated to Rust.


Copying
-------

PolyLog.jl is licenced under the MIT License.
