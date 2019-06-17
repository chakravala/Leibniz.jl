
#   This file is part of Leibniz.jl. It is licensed under the GPL license
#   Leibniz Copyright (C) 2019 Michael Reed

using Grassmann

function (V::Signature{N})(d::Derivation{T,O}) where {N,T,O}
    x=sum([(∂(V,k)^2) for k ∈ 1:N])^div(isodd(O) ? O-1 : O,2)
    isodd(O) ? sum([(x*∂(V,k))*getbasis(V,1<<(k-1)) for k ∈ 1:N]) : x*getbasis(V,0)
end

∂(ω::T) where T<:TensorAlgebra{V} where V = ω⋅V(∇)
d(ω::T) where T<:TensorAlgebra{V} where V = V(∇)∧ω
