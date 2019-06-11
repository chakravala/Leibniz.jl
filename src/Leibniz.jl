module Leibniz

#   This file is part of Leibniz.jl. It is licensed under the GPL license
#   Leibniz Copyright (C) 2019 Michael Reed

using DirectSum, StaticArrays, Requires
using LinearAlgebra, AbstractTensors, AbstractLattices
import Base: *, ^, +, -, /, show
import DirectSum: value

export Differential, Partial, Derivation, d, ∂, ∇, Δ, @operator

abstract type Operator end

struct Differential{D,O,T} <: Operator
    v::T
end

Differential{D}() where D = Differential{D,1}()
Differential{D,O}() where {D,O} = Differential{D,O}(true)
Differential{D,O}(v::T) where {D,O,T} = Differential{D,O,T}(v)

value(d::Differential{D,O,T} where {D,O}) where T = d.v

show(io::IO,d::Differential{D,0,Bool} where D) = print(io,value(d) ? 1 : -1)
show(io::IO,d::Differential{D,0} where D) = print(io,value(d))
show(io::IO,d::Differential{D,1,Bool}) where D = print(io,value(d) ? "" : "-",'∂',DirectSum.subs[D])
show(io::IO,d::Differential{D,1}) where D = print(io,value(d),'∂',DirectSum.subs[D])
show(io::IO,d::Differential{D,O,Bool}) where {D,O} = print(io,value(d) ? "" : "-",'∂',DirectSum.subs[D],DirectSum.sups[O])
show(io::IO,d::Differential{D,O}) where {D,O} = print(io,value(d),'∂',DirectSum.subs[D],DirectSum.sups[O])

+(d::O) where O<:Operator = d
+(r,d::O) where O<:Operator = d+r
+(a::Differential{D,O},b::Differential{D,O}) where {D,O} = (c=a.v+b.v; iszero(c) ? 0 : Differential{D,O}(c))
-(a::Differential{D,O},b::Differential{D,O}) where {D,O} = (c=a.v-b.v; iszero(c) ? 0 : Differential{D,O}(c))
-(d::Differential{D,O}) where {D,O} = Differential{D,O}(-d.v)
*(d::Differential{D,0} where D,r) = r
*(r,d::Differential) = d*r
*(a::Differential{D,1},b::Differential{D,1}) where D = (c=a.v*b.v; iszero(c) ? 0 : Differential{D,2}(c))
*(a::Differential{D,A},b::Differential{D,B}) where {D,A,B} = (c=a.v*b.v; iszero(c) ? 0 : Differential{D,A+B}(c))
*(a::Differential{D,O,Bool},b::I) where {D,O,I<:Number} = isone(b) ? a : Differential{D,O,I}(value(a) ? b : -b)
*(a::Differential{D,O,T},b::I) where {D,O,T,I<:Number} = isone(b) ? a : Differential{D,O}(value(a)*b)
^(d::Differential{D,O},o::T) where {D,O,T<:Integer} = iszero(o) ? 1 : Differential{D,O*o}(value(d)^o)
^(d::Differential{D,O,Bool},o::T) where {D,O,T<:Integer} = iszero(o) ? 1 : Differential{D,O*o}(value(d) ? true : iseven(o))

struct Partial{D,T} <: Operator
    v::T
end

Partial{D}() where D = Partial{D}(true)
Partial{D}(v::T) where {D,T} = Partial{D,T}(v)

value(d::Partial{D,T} where D) where T = d.v

show(io::IO,d::Partial{D,Bool}) where D = print(io,value(d) ? "" : "-",'∂',[DirectSum.subs[k] for k ∈ DirectSum.indices(UInt(D))]...)
show(io::IO,d::Partial{D}) where D = print(io,value(d),'∂',[DirectSum.subs[k] for k ∈ DirectSum.indices(UInt(D))]...)
show(io::IO,d::Partial{0,Bool}) = print(io,value(d) ? 1 : -1)
show(io::IO,d::Partial{0}) = print(io,value(d))
show(io::IO,d::Partial{UInt(0),Bool}) = print(io,value(d) ? 1 : -1)
show(io::IO,d::Partial{UInt(0)}) = print(io,value(d))

indexint(D) = DirectSum.bit2int(DirectSum.indexbits(max(D...),D))

∂(D::T) where T<:Integer = Differential{D}()
∂(D::T...) where T<:Integer = Partial{indexint(D)}()

*(r,d::Partial) = d*r
*(a::Differential{D1,1},b::Differential{D2,1}) where {D1,D2} = (c=a.v*b.v; iszero(c) ? 0 : Partial{indexint((D1,D2))}(c))
*(a::Partial{D,Bool},b::I) where {D,I<:Number} = isone(b) ? a : Partial{D,I}(value(a) ? b : -b)
*(a::Partial{D,T},b::I) where {D,T,I<:Number} = isone(b) ? a : Partial{D}(value(a)*b)
+(a::Partial{D},b::Partial{D}) where D = (c=a.v+b.v; iszero(c) ? 0 : Partial{indexint((D1,D2))}(c))
-(a::Partial{D},b::Partial{D}) where D = (c=a.v-b.v; iszero(c) ? 0 : Partial{indexint((D1,D2))}(c))
#-(d::Partial{D,Bool}) where {D,O} = Partial{D,Bool}(!value(d))
-(d::Partial{D}) where D = Partial{D}(-value(d))

struct OperatorExpr{T} <: Operator
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

+(d::Partial,n::T) where T<:Number = iszero(n) ? d : add(d,n)
+(d::Differential,n::T) where T<:Number = iszero(n) ? d : add(d,n)
+(d::Partial,n::Partial) = add(d,n)
+(d::Differential,n::Partial) = add(d,n)
+(d::Partial,n::Differential) = add(d,n)
+(d::Differential,n::Differential) = add(d,n)
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

*(d::OperatorExpr,n::Differential) = times(d,n)
*(d::OperatorExpr,n::Partial) = times(d,n)
*(d::OperatorExpr,n::OperatorExpr) = OperatorExpr(Expr(:call,:*,d,n))

*(a::Differential,b::Differential{D,O}) where {D,O} = Differential{D,O}(a*OperatorExpr(b.v))
*(a::Differential{D,O},b::OperatorExpr) where {D,O} = Differential{D,O}(value(a)*b.expr)
*(a::Partial,b::Partial{D}) where D = Partial{D}(a*OperatorExpr(b.v))
*(a::Partial{D},b::OperatorExpr) where {D,O} = Partial{D}(value(a)*b.expr)
*(a::Differential,b::Partial{D}) where D = Partial{D}(a*OperatorExpr(b.v))
*(a::Partial,b::Differential{D}) where D = Differential{D}(a*OperatorExpr(b.v))

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

∇ = Derivation(LinearAlgebra.I)
Δ = ∇^2

function __init__()
    @require Grassmann="4df31cd9-4c27-5bea-88d0-e6a7146666d8" include("grassmann.jl")
    @require Reduce="93e0c654-6965-5f22-aba9-9c1ae6b3c259" include("symbolic.jl")
    #@require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" nothing
end

end # module
