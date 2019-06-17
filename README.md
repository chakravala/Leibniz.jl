# Leibniz.jl

*Operator algebras for multivariate differentiable Julia expressions*

[![Build Status](https://travis-ci.org/chakravala/Leibniz.jl.svg?branch=master)](https://travis-ci.org/chakravala/Leibniz.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/xb03dyfvhni6vrj5?svg=true)](https://ci.appveyor.com/project/chakravala/leibniz-jl)
[![Coverage Status](https://coveralls.io/repos/chakravala/Leibniz.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chakravala/Leibniz.jl?branch=master)
[![codecov.io](http://codecov.io/github/chakravala/Leibniz.jl/coverage.svg?branch=master)](http://codecov.io/github/chakravala/Leibniz.jl?branch=master)
[![Gitter](https://badges.gitter.im/Grassmann-jl/community.svg)](https://gitter.im/Grassmann-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
[![Liberapay patrons](https://img.shields.io/liberapay/patrons/chakravala.svg)](https://liberapay.com/chakravala)

Cross-compatibility of [Grassmann.jl](https://github.com/chakravala/Grassmann.jl) with [Reduce.jl](https://github.com/chakravala/Reduce.jl) for multivariable differential operators and tensor field operations.

```Julia
julia> using Reduce, Leibniz, Grassmann
Reduce (Free CSL version, revision 4980), 06-May-19 ...

julia> V = V"3" # load Reduce 1st! otherwise slow
⟨+++⟩

julia> V(∇)
∂₁v₁ + ∂₂v₂ + ∂₃v₃

julia> V(∇^0), V(∇^2)
(1v, (∂₁² + ∂₂² + ∂₃²)v)

julia> V(∇^3)
(∂₁³ + ∂₂²∂₁ + ∂₃²∂₁)v₁ + (∂₁²∂₂ + ∂₂³ + ∂₃²∂₂)v₂ + (∂₁²∂₃ + ∂₂²∂₃ + ∂₃³)v₃

julia> V(∇^4)
((∂₁² + ∂₂² + ∂₃²) ^ 2)v

julia> ∇^2 == Δ
true

julia> ∇, Δ
(∂ₖvₖ, ∂ₖ²v)

```

Generates the tensor algebra of multivariable symmetric Leibniz differentials and interfaces `using Reduce, Grassmann` to provide the `∇,Δ` vector field operators, enabling  mixed-symmetry tensors with arbitrary multivariate `Grassmann` manifolds.

This is an initial undocumented pre-release registration for testing with other packages.
