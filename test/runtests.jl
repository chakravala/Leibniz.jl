using Leibniz, Grassmann
using Test

const V = V"3"

# write your own tests here
@test ∇^2 == Δ
@test value(V(∇)⋅V(∇)) == value((V(∇)^2)(0)(1))
