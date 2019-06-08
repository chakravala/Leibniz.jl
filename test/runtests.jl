using Leibniz,Reduce,Grassmann
using Test

# write your own tests here
@test value(∇(V"3")⋅∇(V"3")).expr == value((∇(V"3")^2)(0)(1)).expr
