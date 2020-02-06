using Leibniz #, Grassmann
using Test

#const V = SubManifold(ℝ^3)

# write your own tests here
@test ∇^2 == Δ
#@test V(∇)⋅V(∇) == V(∇)^2
