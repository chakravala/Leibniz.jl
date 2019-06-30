module Leibniz

#   This file is part of Leibniz.jl. It is licensed under the GPL license
#   Leibniz Copyright (C) 2019 Michael Reed

using DirectSum, StaticArrays #, Requires
using LinearAlgebra, AbstractTensors, AbstractLattices
import Base: *, ^, +, -, /, show, zero
import DirectSum: value, V0, mixedmode, pre

export Differential, Monomial, Derivation, d, ∂, ∇, Δ, @operator

abstract type Operator{V} end #<: TensorAlgebra{V} end

abstract type Polynomial{V,G} <: Operator{V} end

+(d::O) where O<:Operator = d
+(r,d::O) where O<:Operator = d+r

struct Monomial{V,G,D,O,T} <: Polynomial{V,G}
    v::T
end

Monomial{V,G,D}() where {V,G,D} = Monomial{V,G,D}(true)
Monomial{V,G,D}(v::T) where {V,G,D,T} = Monomial{V,G,D,1,T}(v)
Monomial{V,G,D,O}() where {V,G,D,O,T} = Monomial{V,G,D,O}(true)
Monomial{V,G,D,O}(v::T) where {V,G,D,O,T} = Monomial{V,G,D,O,T}(v)

zero(::Monomial) = Monomial{V0,0,0,0}()

value(d::Monomial{V,G,D,T} where {V,G,D}) where T = d.v

sups(O) = O ≠ 1 ? DirectSum.sups[O] : ""

show(io::IO,d::Monomial{V,G,D,O,Bool} where G) where {V,D,O} = print(io,value(d) ? "" : "-",pre[mixedmode(V)>0 ? 4 : 3],[DirectSum.subs[k] for k ∈ DirectSum.shift_indices(V,UInt(D))]...,sups(O))
show(io::IO,d::Monomial{V,G,D,O} where G) where {V,D,O} = print(io,value(d),pre[mixedmode(V)>0 ? 4 : 3],[DirectSum.subs[k] for k ∈ DirectSum.shift_indices(V,UInt(D))]...,sups(O))
show(io::IO,d::Monomial{V,G,D,0,Bool} where {V,G,D}) = print(io,value(d) ? 1 : -1)
show(io::IO,d::Monomial{V,G,0,O,Bool} where {V,G,O}) = print(io,value(d) ? 1 : -1)
show(io::IO,d::Monomial{V,0,D,O,Bool} where {V,D,O}) = print(io,value(d) ? 1 : -1)
show(io::IO,d::Monomial{V,G,D,0} where {V,G,D}) = print(io,value(d))
show(io::IO,d::Monomial{V,G,0} where {V,G}) = print(io,value(d))
show(io::IO,d::Monomial{V,0} where V) = print(io,value(d))
show(io::IO,d::Monomial{V,G,D,UInt(0),Bool} where {V,G,D}) = print(io,value(d) ? 1 : -1)
show(io::IO,d::Monomial{V,G,UInt(0),O,Bool} where {V,G,O}) = print(io,value(d) ? 1 : -1)
show(io::IO,d::Monomial{V,UInt(0),D,O,Bool} where {V,D,O}) = print(io,value(d) ? 1 : -1)
show(io::IO,d::Monomial{V,G,D,UInt(0)} where {V,G,D}) = print(io,value(d))
show(io::IO,d::Monomial{V,G,UInt(0)} where {V,G}) = print(io,value(d))
show(io::IO,d::Monomial{V,UInt(0)} where V) = print(io,value(d))

indexint(D) = DirectSum.bit2int(DirectSum.indexbits(max(D...),D))

∂(D::T...) where T<:Integer = Monomial{V0,length(D),indexint(D)}()
∂(V::S,D::T...) where {S<:VectorBundle,T<:Integer} = Monomial{V,length(D),indexint(D)}()

*(r,d::Monomial) = d*r
*(d::Monomial{V,G,D,0} where {V,G,D},r) = r
*(d::Monomial{V,G,0,O} where {V,G,O},r) = r
*(d::Monomial{V,0,D,O} where {V,D,O},r) = r
*(a::Monomial{V,1,D,O1},b::Monomial{V,1,D,O2}) where {V,D,O1,O2} = (c=a.v*b.v; iszero(c) ? 0 : Monomial{V,2,D,O1+O2}(c))
*(a::Monomial{V,1,D,1},b::Monomial{V,1,D,1}) where {V,D,O1,O2} = (c=a.v*b.v; iszero(c) ? 0 : Monomial{V,2,D,2}(c))
*(a::Monomial{V,1,D1,1},b::Monomial{V,1,D2,1}) where {V,D1,D2} = (c=a.v*b.v; iszero(c) ? 0 : Monomial{V,2,D1|D2}(c))
*(a::Monomial{V,G,D,O,Bool},b::I) where {V,G,D,O,I<:Number} = isone(b) ? a : Monomial{V,G,D,O,I}(value(a) ? b : -b)
*(a::Monomial{V,G,D,O,T},b::I) where {V,G,D,O,T,I<:Number} = isone(b) ? a : Monomial{V,G,D,O}(value(a)*b)
+(a::Monomial{V,G,D,O},b::Monomial{V,G,D,O}) where {V,G,D,O} = (c=a.v+b.v; iszero(c) ? 0 : Monomial{V,G,D,O}(c))
-(a::Monomial{V,G,D,O},b::Monomial{V,G,D,O}) where {V,G,D,O} = (c=a.v-b.v; iszero(c) ? 0 : Monomial{V,G,D,O}(c))
#-(d::Monomial{V,G,D,O,Bool}) where {V,G,D,O} = Monomial{V,G,D,O,Bool}(!value(d))
-(d::Monomial{V,G,D,O}) where {V,G,D,O} = Monomial{V,G,D,O}(-value(d))
^(d::Monomial{V,G,D,O},o::T) where {V,G,D,O,T<:Integer} = iszero(o) ? 1 : Monomial{V,G+(O*o),D,O*o}(value(d)^o)
^(d::Monomial{V,G,D,O,Bool},o::T) where {V,G,D,O,T<:Integer} = iszero(o) ? 1 : Monomial{V,G+(O*o),D,O*o}(value(d) ? true : iseven(o))

struct OperatorExpr{T} <: Operator{T}
    expr::T
end

macro operator(expr)
    OperatorExpr(expr)
end

show(io::IO,d::OperatorExpr) = print(io,'(',d.expr,')')

add(d,n) = OperatorExpr(Expr(:call,:+,d,n))

function plus(d::OperatorExpr{T},n) where T
    iszero(n) && (return d)
    if T == Expr
        if d.expr.head == :call
            if d.expr.args[1] == :+
                return OperatorExpr(Expr(:call,:+,push!(copy(d.expr.args[2:end]),n)...))
            elseif d.expr.args[1] == :- && length(d.expr.args) == 2 && d.expr.args[2] == n
                return 0
            else
                return OperatorExpr(Expr(:call,:+,d.expr,n))
            end
        else
            throw(error("Operator expression not implemented"))
        end
    else
        OperatorExpr(d.expr+n)
    end
end

+(d::Monomial,n::T) where T<:Number = iszero(n) ? d : add(d,n)
+(d::Monomial,n::Monomial) = add(d,n)
+(d::OperatorExpr,n) = plus(d,n)
+(d::OperatorExpr,n::O) where O<:Operator = plus(d,n)
-(d::OperatorExpr) = OperatorExpr(Expr(:call,:-,d))

#add(d,n) = OperatorExpr(Expr(:call,:+,d,n))

function times(d::OperatorExpr{T},n) where T
    iszero(n) && (return 0)
    isone(n) && (return d)
    if T == Expr
        if d.expr.head == :call
            if d.expr.args[1] ∈ (:+,:-)
                return OperatorExpr(Expr(:call,:+,(d.expr.args[2:end] .* Ref(n))...))
            elseif d.expr.args[1] == :*
                (d.expr.args[2]*n)*d.expr.args[3] + d.expr.args[2]*(d.expr.args[3]*n)
            elseif d.expr.args[1] == :/
                (d.expr.args[2]*n)/d.expr.args[3] - (d.expr.args[2]*(d.expr.args[3]*n))/(d.expr.args[3]^2)
            else
                return OperatorExpr(Expr(:call,:*,d.expr,n))
            end
        else
            throw(error("Operator expression not implemented"))
        end
    else
        OperatorExpr(d.expr*n)
    end
end

*(d::OperatorExpr,n::Monomial) = times(d,n)
*(d::OperatorExpr,n::OperatorExpr) = OperatorExpr(Expr(:call,:*,d,n))

*(a::Monomial{V,1},b::Monomial{V,1,D,O}) where {V,D,O} = Monomial{V,1,D,O}(a*OperatorExpr(b.v))
*(a::Monomial,b::Monomial{V,G,D,O}) where {V,G,D,O} = Monomial{V,G,D,O}(a*OperatorExpr(b.v))
*(a::Monomial{V,G,D,O},b::OperatorExpr) where {V,G,D,O} = Monomial{V,G,D,O}(value(a)*b.expr)

^(d::OperatorExpr,n::T) where T<:Integer = iszero(n) ? 1 : isone(n) ? d : OperatorExpr(Expr(:call,:^,d,n))

## generic

Base.signbit(::O) where O<:Operator = false
Base.abs(d::O) where O<:Operator = d

struct Derivation{T,O}
    v::UniformScaling{T}
end

Derivation{T}(v::UniformScaling{T}) where T = Derivation{T,1}(v)
Derivation(v::UniformScaling{T}) where T = Derivation{T}(v)

show(io::IO,v::Derivation{Bool,O}) where O = print(io,(v.v.λ ? "" : "-"),"∂ₖ",O==1 ? "" : DirectSum.sups[O],"v",isodd(O) ? "ₖ" : "")
show(io::IO,v::Derivation{T,O}) where {T,O} = print(io,v.v.λ,"∂ₖ",O==1 ? "" : DirectSum.sups[O],"v",isodd(O) ? "ₖ" : "")

-(v::Derivation{Bool,O}) where {T,O} = Derivation{Bool,O}(UniformScaling{Bool}(!v.v.λ))
-(v::Derivation{T,O}) where {T,O} = Derivation{T,O}(UniformScaling{T}(-v.v.λ))

function ^(v::Derivation{T,O},n::S) where {T,O,S<:Integer}
    x = T<:Bool ? (isodd(n) ? v.v.λ : true ) : v.v.λ^n
    t = typeof(x)
    Derivation{t,O*n}(UniformScaling{t}(x))
end

for op ∈ (:+,:-,:*)
    @eval begin
        $op(a::Derivation{A,O},b::Derivation{B,O}) where {A,B,O} = Derivation{promote_type(A,B),O}($op(a.v,b.v))
        $op(a::Derivation{A,O},b::B) where {A,B<:Number,O} = Derivation{promote_type(A,B),O}($op(a.v,b))
        $op(a::A,b::Derivation{B,O}) where {A<:Number,B,O} = Derivation{promote_type(A,B),O}($op(a,b.v))
    end
end

unitype(::UniformScaling{T}) where T = T

/(a::Derivation{A,O},b::Derivation{B,O}) where {A,B,O} = (x=a.v/b.v; Derivation{unitype(x),O}(x))
/(a::Derivation{A,O},b::B) where {A,B<:Number,O} = (x=a.v/b; Derivation{unitype(x),O}(x))
#/(a::A,b::Derivation{B,O}) where {A<:Number,B,O} = (x=a/b.v; Derivation{typeof(x),O}(x))

import AbstractLattices: ∧, ∨
import LinearAlgebra: dot, cross

for op ∈ (:+,:-,:*,:/,:∧,:∨,:dot,:cross)
    @eval begin
        $op(a::Derivation,b::B) where B<:TensorAlgebra{V} where V = $op(V(a),b)
        $op(a::A,b::Derivation) where A<:TensorAlgebra{V} where V = $op(a,V(b))
    end
end

const ∇ = Derivation(LinearAlgebra.I)
const Δ = ∇^2

function d end

#=function __init__()
    @require Grassmann="4df31cd9-4c27-5bea-88d0-e6a7146666d8" include("grassmann.jl")
    @require Reduce="93e0c654-6965-5f22-aba9-9c1ae6b3c259" include("symbolic.jl")
    #@require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" nothing
end=#

end # module
