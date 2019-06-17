using Leibniz # Grassmann
using Test

# write your own tests here
@test ∇^2 == Δ
#@test value(V(∇)⋅V(∇)).expr == value((V(∇)^2)(0)(1)).expr
