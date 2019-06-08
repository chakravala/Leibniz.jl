module Leibniz

#   This file is part of Leibniz.jl. It is licensed under the GPL license
#   Leibniz Copyright (C) 2019 Michael Reed

using DirectSum, StaticArrays, Requires
import Base: *, ^, +, -, /, show

export Differential, Partial, ∂

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

struct OperatorExpr <: Operator
    expr
end

show(io::IO,d::OperatorExpr) = print(io,'(',d.expr,')')

add(d,n) = OperatorExpr(Expr(:call,:+,d,n))

function plus(d::OperatorExpr,n)
    iszero(n) && (return d)
    if typeof(d.expr) == Expr
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
        OperatorExpr(d+n)
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

## generic

Base.signbit(::O) where O<:Operator = false
Base.abs(d::O) where O<:Operator = d

function __init__()
    @require Grassmann="4df31cd9-4c27-5bea-88d0-e6a7146666d8" include("grassmann.jl")
    @require Reduce="93e0c654-6965-5f22-aba9-9c1ae6b3c259" include("symbolic.jl")
    #@require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" nothing
end

end # module
