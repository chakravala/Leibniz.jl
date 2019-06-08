
#   This file is part of Leibniz.jl. It is licensed under the GPL license
#   Leibniz Copyright (C) 2019 Michael Reed

using Reduce

Reduce.RExpr(::Differential{D,0}) where D = RExpr("1")
Reduce.RExpr(::Differential{D,1}) where D = RExpr("d_$D")
Reduce.RExpr(::Differential{D,O}) where {D,O} = RExpr("d_$D^$O")

for RE ∈ (RExpr,Expr,Symbol)
    @eval begin
        *(d::Differential{D,1},r::$RE) where D = Algebra.df(r,RExpr("v$D"))
        *(d::Differential{D,O},r::$RE) where {D,O} = Algebra.df(r,RExpr("v$D"),RExpr("$O"))
    end
end
