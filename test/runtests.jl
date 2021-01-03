using Leibniz, Grassmann
using Test

# write your own tests here
@test ∇^2 == Δ
@test ℝ3(∇)⋅ℝ3(∇) == ℝ3(∇)^2
