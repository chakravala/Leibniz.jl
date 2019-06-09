module Leibniz

#   This file is part of Leibniz.jl. It is licensed under the GPL license
#   Leibniz Copyright (C) 2019 Michael Reed

using DirectSum, StaticArrays, Requires
using LinearAlgebra, AbstractTensors, AbstractLattices
import Base: *, ^, +, -, /, show

export Differential, Partial, Derivation, d, ∂, ∇, Δ

abstract type Operator end

struct Differential{D,O} <: Operator end
Differential{D}() where D = Differential{D,1}()

show(io::IO,::Differential{D,0} where D) = print(io,1)
show(io::IO,::Differential{D,1}) where D = print(io,'∂',DirectSum.subs[D])
show(io::IO,::Differential{D,O}) where {D,O} = print(io,'∂',DirectSum.subs[D],DirectSum.sups[O])

+(d::Differential) = d
+(r,d::O) where O<:Operator = d+r
*(d::Differential{D,0} where D,r) = r
*(r,d::Differential) = d*r
*(::Differential{D,1},::Differential{D,1}) where D = Differential{D,2}()
*(::Differential{D,O1},::Differential{D,O2}) where {D,O1,O2} = Differential{D,O1+O2}()
^(::Differential{D,O},o::T) where {D,O,T<:Integer} = Differential{D,O*o}()

struct Partial{D} <: Operator end

show(io::IO,::Partial{D}) where D = print(io,'∂',[DirectSum.subs[k] for k ∈ DirectSum.indices(UInt(D))]...)
show(io::IO,::Partial{0}) = print(io,1)
show(io::IO,::Partial{UInt(0)}) = print(io,1)

indexint(D) = DirectSum.bit2int(DirectSum.indexbits(max(D...),D))

∂(D::T) where T<:Integer = Differential{D}()
∂(D::T...) where T<:Integer = Partial{indexint(D)}()

*(::Differential{D1,1},::Differential{D2,1}) where {D1,D2} = ∂(D1,D2)

struct OperatorExpr{T} <: Operator
    expr::T
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
-(d::O) where O<:Operator = OperatorExpr(Expr(:call,:-,d))

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
*(n::Differential,d::OperatorExpr) = times(d,n)
*(n::Partial,d::OperatorExpr) = times(d,n)

## generic

Base.signbit(::O) where O<:Operator = false
Base.abs(d::O) where O<:Operator = d

struct Derivation{T,O}
    v::UniformScaling{T}
end

Derivation{T}(v::UniformScaling{T}) where T = Derivation{T,1}(v)
Derivation(v::UniformScaling{T}) where T = Derivation{T}(v)

show(io::IO,v::Derivation{Bool,O}) where O = print(io,(v.v.λ ? "" : "-"),"∂ₖ",O==1 ? "" : DirectSum.sups[O],"vₖ")
show(io::IO,v::Derivation{T,O}) where {T,O} = print(io,v.v.λ,"∂ₖ",O==1 ? "" : DirectSum.sups[O],"vₖ")

-(v::Derivation{Bool,O}) where {T,O} = Derivation{Bool,O}(LinearAlgebra.UniformScaling{Bool}(!v.v.λ))
-(v::Derivation{T,O}) where {T,O} = Derivation{T,O}(LinearAlgebra.UniformScaling{T}(-v.v.λ))

function ^(v::Derivation{T,O},n::S) where {T,O,S<:Integer}
    x = T<:Bool ? (isodd(n) ? v.v.λ : true ) : v.v.λ^n
    t = typeof(x)
    Derivation{t,O*n}(LinearAlgebra.UniformScaling{t}(x))
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

include("grassmann.jl")

function __init__()
    #@require Grassmann="4df31cd9-4c27-5bea-88d0-e6a7146666d8" include("grassmann.jl")
    @require Reduce="93e0c654-6965-5f22-aba9-9c1ae6b3c259" include("symbolic.jl")
    #@require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" nothing
end

end # module
