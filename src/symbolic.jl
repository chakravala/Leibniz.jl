
#   This file is part of Leibniz.jl. It is licensed under the GPL license
#   Leibniz Copyright (C) 2019 Michael Reed

using Reduce

Reduce.RExpr(::Monomial{V,G,D,0} where {V,G}) where D = RExpr("1")
Reduce.RExpr(::Monomial{V,G,D,1} where {V,G}) where D = RExpr("d_$D")
Reduce.RExpr(::Monomial{V,G,D,O} where {V,G}) where {D,O} = RExpr("d_$D^$O")

for RE ∈ (RExpr,Expr,Symbol)
    @eval begin
        *(d::Monomial{V,1,D,1} where V,r::$RE) where D = Algebra.df(r,RExpr("v$D"))
        *(d::Monomial{V,O,D,O} where V,r::$RE) where {D,O} = Algebra.df(r,RExpr("v$D"),RExpr("$O"))
    end
end
