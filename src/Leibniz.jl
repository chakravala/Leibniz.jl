module Leibniz

#   This file is part of Leibniz.jl. It is licensed under the AGPL license
#   Leibniz Copyright (C) 2019 Michael Reed

using LinearAlgebra, AbstractTensors
export Manifold, Differential, Derivation, d, ∂, ∇, Δ
import Base: getindex, convert, @pure, +, *, ∪, ∩, ⊆, ⊇, ==, show, zero
import LinearAlgebra: det, rank

## Manifold{N}

import AbstractTensors: TensorAlgebra, Manifold, TensorGraded, scalar, isscalar, involute
import AbstractTensors: vector, isvector, bivector, isbivector, volume, isvolume, ⋆, mdims
import AbstractTensors: value, valuetype, interop, interform, even, odd, isnull, norm
import AbstractTensors: TupleVector, Values, Variables, FixedVector, SVector, MVector
import AbstractTensors: basis, complementleft, complementlefthodge, unit, involute, clifford
abstract type TensorTerm{V,G} <: TensorGraded{V,G} end

## utilities

include("utilities.jl")

#="""
    floatprecision(s)

Set float precision for display Float64 coefficents.

Float coefficients `f` are printed as `@sprintf(s,f)`.

If `s == ""` (default), then `@sprintf` is not used.
"""
const floatprecision = ( () -> begin
        gs::String = ""
        return (tf=gs)->(gs≠tf && (gs=tf); return gs)
    end)()
export floatprecision

macro fprintf()
    s = floatprecision()
    isempty(s) ? :(m.v) : :(Printf.@sprintf($s,m.v))
end=#

# symbolic print types

parval = (Expr,Complex,Rational,TensorAlgebra)

# number fields

const Fields = (Real,Complex)
const Field = Fields[1]
const ExprField = Union{Expr,Symbol}

extend_field(Field=Field) = (global parval = (parval...,Field))

for T ∈ Fields
    @eval begin
        Base.:(==)(a::T,b::TensorTerm{V,G} where V) where {T<:$T,G} = G==0 ? a==value(b) : 0==a==value(b)
        Base.:(==)(a::TensorTerm{V,G} where V,b::T) where {T<:$T,G} = G==0 ? value(a)==b : 0==value(a)==b
    end
end

Base.:(==)(a::TensorTerm,b::TensorTerm) = 0 == value(a) == value(b)

for T ∈ (Fields...,Symbol,Expr)
    @eval begin
        Base.isapprox(a::S,b::T) where {S<:TensorAlgebra,T<:$T} = Base.isapprox(a,Simplex{Manifold(a)}(b))
        Base.isapprox(a::S,b::T) where {S<:$T,T<:TensorAlgebra} = Base.isapprox(b,a)
    end
end

## fundamentals

"""
    getbasis(V::Manifold,v)

Fetch a specific `SubManifold{G,V}` element from an optimal `SubAlgebra{V}` selection.
"""
@inline getbasis(V,b) = getbasis(V,UInt(b))

Base.one(V::T) where T<:TensorGraded = one(Manifold(V))
Base.zero(V::T) where T<:TensorGraded = zero(Manifold(V))

@pure g_one(::Type{T}) where T = one(T)
@pure g_zero(::Type{T}) where T = zero(T)

## Derivation

struct Derivation{T,O}
    v::UniformScaling{T}
end

Derivation{T}(v::UniformScaling{T}) where T = Derivation{T,1}(v)
Derivation(v::UniformScaling{T}) where T = Derivation{T}(v)

show(io::IO,v::Derivation{Bool,O}) where O = print(io,(v.v.λ ? "" : "-"),"∂ₖ",O==1 ? "" : AbstractTensors.sups[O],"v",isodd(O) ? "ₖ" : "")
show(io::IO,v::Derivation{T,O}) where {T,O} = print(io,v.v.λ,"∂ₖ",O==1 ? "" : AbstractTensors.sups[O],"v",isodd(O) ? "ₖ" : "")

Base.:-(v::Derivation{Bool,O}) where {T,O} = Derivation{Bool,O}(UniformScaling{Bool}(!v.v.λ))
Base.:-(v::Derivation{T,O}) where {T,O} = Derivation{T,O}(UniformScaling{T}(-v.v.λ))

function Base.:^(v::Derivation{T,O},n::S) where {T,O,S<:Integer}
    x = T<:Bool ? (isodd(n) ? v.v.λ : true ) : v.v.λ^n
    t = typeof(x)
    Derivation{t,O*n}(UniformScaling{t}(x))
end

for op ∈ (:(Base.:+),:(Base.:-),:(Base.:*))
    @eval begin
        $op(a::Derivation{A,O},b::Derivation{B,O}) where {A,B,O} = Derivation{promote_type(A,B),O}($op(a.v,b.v))
        $op(a::Derivation{A,O},b::B) where {A,B<:Number,O} = Derivation{promote_type(A,B),O}($op(a.v,b))
        $op(a::A,b::Derivation{B,O}) where {A<:Number,B,O} = Derivation{promote_type(A,B),O}($op(a,b.v))
    end
end

unitype(::UniformScaling{T}) where T = T

Base.:/(a::Derivation{A,O},b::Derivation{B,O}) where {A,B,O} = (x=Base.:/(a.v,b.v); Derivation{unitype(x),O}(x))
Base.:/(a::Derivation{A,O},b::B) where {A,B<:Number,O} = (x=Base.:/(a.v,b); Derivation{unitype(x),O}(x))
#Base.:/(a::A,b::Derivation{B,O}) where {A<:Number,B,O} = (x=Base.:/(a,b.v); Derivation{typeof(x),O}(x))
Base.:\(a::Derivation{A,O},b::Derivation{B,O}) where {A,B,O} = (x=a.v\b.v; Derivation{unitype(x),O}(x))
Base.:\(a::A,b::Derivation{B,O}) where {A<:Number,B,O} = (x=a\b.v; Derivation{unitype(x),O}(x))

import AbstractTensors: ∧, ∨
import LinearAlgebra: dot, cross

for op ∈ (:(Base.:+),:(Base.:-),:(Base.:*),:(Base.:/),:(Base.:\),:∧,:∨,:dot,:cross)
    @eval begin
        $op(a::Derivation,b::B) where B<:TensorAlgebra = $op(Manifold(b)(a),b)
        $op(a::A,b::Derivation) where A<:TensorAlgebra = $op(a,Manifold(a)(b))
    end
end

const ∇ = Derivation(LinearAlgebra.I)
const Δ = ∇^2

function d end
function ∂ end

include("generic.jl")
include("operations.jl")
include("indices.jl")

bladeindex(cache_limit,one(UInt))
basisindex(cache_limit,one(UInt))

indexbasis(Int((sparse_limit+cache_limit)/2),1)

#=function __init__()
    @require Reduce="93e0c654-6965-5f22-aba9-9c1ae6b3c259" include("symbolic.jl")
    #@require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" nothing
end=#

end # module
