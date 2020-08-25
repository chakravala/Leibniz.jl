
#   This file is part of Leibniz.jl. It is licensed under the AGPL license
#   Leibniz Copyright (C) 2019 Michael Reed

export bits, basis, grade, order, options, metric, polymode, dyadmode, diffmode, diffvars
export valuetype, value, hasinf, hasorigin, isorigin, norm, indices, tangent, isbasis, ≅

@pure grade(V::M) where M<:Manifold{N} where N = N-(isdyadic(V) ? 2 : 1)*diffvars(V)
@pure grade(m::T) where T<:Real = 0
@pure order(m) = 0
@pure order(V::M) where M<:Manifold = diffvars(V)
@pure options(::Int) = 0
@pure metric(::Int) = zero(UInt)
@pure metric(V::M,b::UInt) where M<:Manifold = PROD(V[indices(b)])
@pure polymode(::Int) = true
@pure dyadmode(::Int) = 0
@pure diffmode(::Int) = 0
@pure diffvars(::Int) = 0
for mode ∈ (:options,:polymode,:dyadmode,:diffmode,:diffvars)
    @eval @pure $mode(t::T) where T<:TensorAlgebra = $mode(Manifold(t))
end

@pure ≅(a,b) = grade(a) == grade(b) && order(a) == order(b) && diffmode(a) == diffmode(b)

export isdyadic, isdual, istangent
const mixedmode = dyadmode
@pure isdyadic(t::T) where T<:TensorAlgebra = dyadmode(Manifold(t))<0
@pure isdual(t::T) where T<:TensorAlgebra = dyadmode(Manifold(t))>0
@pure istangent(t::T) where T<:TensorAlgebra = diffvars(Manifold(t))≠0

@pure isbasis(x) = false
@pure isdyadic(::Int) = false
@pure isdual(::Int) = false
@pure istangent(::Int) = false

@inline value_diff(m::T) where T<:TensorTerm = (v=value(m);istensor(v) ? v : m)

for T ∈ (:T,:(Type{T}))
    @eval @pure isbasis(::$T) where T<:TensorAlgebra = false
end
@pure UInt(m::T) where T<:TensorTerm = UInt(basis(m))

@pure hasconformal(V) = hasinf(V) && hasorigin(V)
@pure hasinf(::Int) = false
@pure hasinf(t::M) where M<:Manifold = hasinf(Manifold(t))
#@pure hasinf(m::T) where T<:TensorAlgebra = hasinf(Manifold(m))
@pure hasorigin(::Int) = false
@pure hasorigin(t::M) where M<:Manifold = hasorigin(Manifold(t))
#@pure hasorigin(m::T) where T<:TensorAlgebra = hasorigin(Manifold(m))

@pure hasorigin(V,B::UInt) = hasinf(V) ? (Bits(2)&B)==Bits(2) : isodd(B)

@pure function hasinf(V,A::UInt,B::UInt)
    hasconformal(V) && (isodd(A) || isodd(B))
end
@pure function hasorigin(V,A::UInt,B::UInt)
    hasconformal(V) && (hasorigin(V,A) || hasorigin(V,B))
end

@pure function hasinf2(V,A::UInt,B::UInt)
    hasconformal(V) && isodd(A) && isodd(B)
end
@pure function hasorigin2(V,A::UInt,B::UInt)
    hasconformal(V) && hasorigin(V,A) && hasorigin(V,B)
end

@pure function hasorigininf(V,A::UInt,B::UInt)
    hasconformal(V) && hasorigin(V,A) && isodd(B) && !hasorigin(V,B) && !isodd(A)
end
@pure function hasinforigin(V,A::UInt,B::UInt)
    hasconformal(V) && isodd(A) && hasorigin(V,B) && !isodd(B) && !hasorigin(V,A)
end

@pure function hasinf2origin(V,A::UInt,B::UInt)
    hasinf2(V,A,B) && hasorigin(V,A,B)
end
@pure function hasorigin2inf(V,A::UInt,B::UInt)
    hasorigin2(V,A,B) && hasinf(V,A,B)
end

@pure diffmask(::Int) = zero(UInt)
@pure function diffmask(V::M) where M<:Manifold
    d = diffvars(V)
    if isdyadic(V)
        v = ((one(UInt)<<d)-1)<<(mdims(V)-2d)
        w = ((one(UInt)<<d)-1)<<(mdims(V)-d)
        return d<0 ? (typemax(UInt)-v,typemax(UInt)-w) : (v,w)
    end
    v = ((one(UInt)<<d)-1)<<(mdims(V)-d)
    d<0 ? typemax(UInt)-v : v
end

@pure function symmetricsplit(V,a)
    sm,dm = symmetricmask(V,a),diffmask(V)
    isdyadic(V) ? (sm&dm[1],sm&dm[2]) : sm
end

@pure function symmetricmask(V,a)
    d = diffmask(V)
    a&(isdyadic(V) ? |(d...) : d)
end

@pure function symmetricmask(V,a,b)
    d = diffmask(V)
    D = isdyadic(V) ? |(d...) : d
    aD,bD = (a&D),(b&D)
    return a&~D, b&~D, aD|bD, aD&bD
end

@pure function diffcheck(V,A::UInt,B::UInt)
    d,db = diffvars(V),diffmask(V)
    v = isdyadic(V) ? db[1]|db[2] : db
    hi = hasinf2(V,A,B) && !hasorigin(V,A,B)
    ho = hasorigin2(V,A,B) && !hasinf(V,A,B)
    (hi || ho) || (d≠0 && count_ones(A&v)+count_ones(B&v)>diffmode(V))
end

## functors

@pure loworder(N::Int) = N
function supermanifold end

## adjoint parities

@pure parityreverse(G) = isodd(Int((G-1)*G/2))
@pure parityinvolute(G) = isodd(G)
@pure parityclifford(G) = parityreverse(G)⊻parityinvolute(G)
const parityconj = parityreverse

## reverse

import Base: reverse, ~
export involute, clifford

@pure grade_basis(V::Int,B) = B&(one(UInt)<<V-1)
@pure grade_basis(V,B) = B&(one(UInt)<<grade(V)-1)
@pure grade(V,B) = count_ones(grade_basis(V,B))

@doc """
    ~(ω::TensorAlgebra)

Reverse of an element: ~ω = (-1)^(grade(ω)*(grade(ω)-1)/2)*ω
""" Base.conj
#reverse(a::UniformScaling{Bool}) = UniformScaling(!a.λ)
#reverse(a::UniformScaling{T}) where T<:Field = UniformScaling(-a.λ)

"""
    reverse(ω::TensorAlgebra)

Reverse of an element: ~ω = (-1)^(grade(ω)*(grade(ω)-1)/2)*ω
"""
@inline Base.:~(b::TensorAlgebra) = Base.conj(b)
#@inline ~(b::UniformScaling) = reverse(b)

@doc """
    involute(ω::TensorAlgebra)

Involute of an element: ~ω = (-1)^grade(ω)*ω
""" involute

@doc """
    clifford(ω::TensorAlgebra)

Clifford conjugate of an element: clifford(ω) = involute(conj(ω))
""" clifford

odd(t::T) where T<:TensorGraded{V,G} where {V,G} = parityinvolute(G) ? t : zero(V)
even(t::T) where T<:TensorGraded{V,G} where {V,G} = parityinvolute(G) ? zero(V) : t

"""
    imag(ω::TensorAlgebra)

The `imag` part `(ω-(~ω))/2` is defined by `abs2(imag(ω)) == -(imag(ω)^2)`.
"""
Base.imag(t::T) where T<:TensorGraded{V,G} where {V,G} = parityreverse(G) ? t : zero(V)

"""
real(ω::TensorAlgebra)

The `real` part `(ω+(~ω))/2` is defined by `abs2(real(ω)) == real(ω)^2`.
"""
Base.real(t::T) where T<:TensorGraded{V,G} where {V,G} = parityreverse(G) ? zero(V) : t

Base.isfinite(b::T) where T<:TensorTerm = isfinite(value(b))

# comparison (special case for scalars)

Base.isless(a::T,b::S) where {T<:TensorTerm{V,0},S<:TensorTerm{W,0}} where {V,W} = isless(value(a),value(b))
Base.isless(a::T,b) where T<:TensorTerm{V,0} where V = isless(value(a),b)
Base.isless(a,b::T) where T<:TensorTerm{V,0} where V = isless(a,value(b))
Base.:<=(x::T,y::S) where {T<:TensorTerm{V,0},S<:TensorTerm{W,0}} where {V,W} = isless(x,y) | (x == y)
Base.:<=(x::T,y) where T<:TensorTerm{V,0} where V = isless(x,y) | (x == y)
Base.:<=(x,y::T) where T<:TensorTerm{V,0} where V = isless(x,y) | (x == y)
