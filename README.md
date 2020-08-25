# Leibniz.jl

*Bit entanglements for tensor algebra derivations and hypergraphs*

[![Build Status](https://travis-ci.org/chakravala/Leibniz.jl.svg?branch=master)](https://travis-ci.org/chakravala/Leibniz.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/xb03dyfvhni6vrj5?svg=true)](https://ci.appveyor.com/project/chakravala/leibniz-jl)
[![Coverage Status](https://coveralls.io/repos/chakravala/Leibniz.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chakravala/Leibniz.jl?branch=master)
[![codecov.io](http://codecov.io/github/chakravala/Leibniz.jl/coverage.svg?branch=master)](http://codecov.io/github/chakravala/Leibniz.jl?branch=master)
[![Gitter](https://badges.gitter.im/Grassmann-jl/community.svg)](https://gitter.im/Grassmann-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
[![Liberapay patrons](https://img.shields.io/liberapay/patrons/chakravala.svg)](https://liberapay.com/chakravala)

Although intended for compatibility use with the [Grassmann.jl](https://github.com/chakravala/Grassmann.jl) package for multivariable differential operators and tensor field operations, `Leibniz` can be used independently.

### Extended dual index printing with full alphanumeric characters #62'

To help provide a commonly shared and readable indexing to the user, some print methods are provided:
```julia
julia> Leibniz.printindices(stdout,Leibniz.indices(UInt(2^62-1)),false,"v")
v₁₂₃₄₅₆₇₈₉₀abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ

julia> Leibniz.printindices(stdout,Leibniz.indices(UInt(2^62-1)),false,"w")
w¹²³⁴⁵⁶⁷⁸⁹⁰ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz
```
An application of this is in `Grassmann` and `DirectSum`, where dual indexing is used.

# Derivation

Generates the tensor algebra of multivariable symmetric Leibniz differentials and interfaces `using Reduce, Grassmann` to provide the `∇,Δ` vector field operators, enabling  mixed-symmetry tensors with arbitrary multivariate `Grassmann` manifolds.

```Julia
julia> using Leibniz, Grassmann

julia> V = tangent(ℝ^3,4,3)
⟨+++⟩

julia> V(∇)
∂₁v₁ + ∂₂v₂ + ∂₃v₃

julia> V(∇^2)
0 + 1∂₁∂₁ + 1∂₂∂₂ + 1∂₃∂₃

julia> V(∇^3)
0 + 1∂₁∂₁∂₁v₁ + 1∂₂∂₂∂₂v₂ + 1∂₃∂₃∂₃v₃ + 1∂₂∂₁₂v₁ + 1∂₃∂₁₃v₁ + 1∂₁∂₁₂v₂ + 1∂₃∂₂₃v₂ + 1∂₁∂₁₃v₃ + 1∂₂∂₂₃v₃

julia> V(∇^4)
0.0 + 1∂₁∂₁∂₁∂₁ + 1∂₂∂₂∂₂∂₂ + 1∂₃∂₃∂₃∂₃ + 2∂₁₂∂₁₂ + 2∂₁₃∂₁₃ + 2∂₂₃∂₂₃

julia> ∇^2 == Δ
true

julia> ∇, Δ
(∂ₖvₖ, ∂ₖ²v)
```

This is an initial undocumented pre-release registration for testing with other packages.
