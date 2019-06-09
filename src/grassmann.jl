
#   This file is part of Leibniz.jl. It is licensed under the GPL license
#   Leibniz Copyright (C) 2019 Michael Reed

using Grassmann

function (V::Signature)(d::Derivation{T,O}) where {T,O,S}
    sum([(∂(k)^O)*getbasis(V,1<<(k-1)) for k ∈ 1:ndims(V)])
end

∂(ω::T) where T<:TensorAlgebra{V} where V = ω⋅V(∇)
d(ω::T) where T<:TensorAlgebra{V} where V = V(∇)∧ω
