
#   This file is part of Leibniz.jl. It is licensed under the AGPL license
#   Leibniz Copyright (C) 2019 Michael Reed

export basis, grade, order, options, metric, polymode, dyadmode, diffmode, diffvars
export valuetype, value, hasinf, hasorigin, isorigin, norm, indices, tangent, isbasis, ≅
export pseudograde

@pure grade(::Type{<:TensorGraded{V,G}}) where {V,G} = G-(isdyadic(V) ? 2 : 1)*diffvars(V)
@pure pseudograde(::Type{<:TensorGraded{V,G}}) where {V,G} = mdims(V)-G-(isdyadic(V) ? 2 : 1)*diffvars(V)
@pure pseudograde(V::M) where M<:Manifold = mdims(V)-rank(V)-(isdyadic(V) ? 2 : 1)*diffvars(V)
@pure grade(V::M) where M<:Manifold = rank(V)-(isdyadic(V) ? 2 : 1)*diffvars(V)
@pure grade(m::T) where T<:Real = 0
@pure order(m) = 0
@pure order(V::M) where M<:Manifold = diffvars(V)
@pure options(::Int) = 0
@pure metric(::Int) = zero(UInt)
@pure metric(V::M,b::UInt) where M<:Manifold = PROD(V[indices(b)])
@pure polymode(::Int) = true
@pure dyadmode(::Int) = 0
@pure diffmode(::Int) = 0
@pure diffvars(::Int) = 0
for mode ∈ (:options,:polymode,:dyadmode,:diffmode,:diffvars)
    @eval begin
        @pure $mode(t::T) where T<:TensorAlgebra = $mode(Manifold(t))
        @pure $mode(t::Type{<:TensorAlgebra}) = $mode(Manifold(t))
    end
end

@pure ≅(a,b) = grade(a) == grade(b) && order(a) == order(b) && diffmode(a) == diffmode(b)

export isdyadic, isdual, istangent
const mixedmode,antigrade = dyadmode,pseudograde
@pure isdyadic(t::Type{<:TensorAlgebra}) = dyadmode(Manifold(t))<0
@pure isdual(t::Type{<:TensorAlgebra}) = dyadmode(Manifold(t))>0
@pure istangent(t::Type{<:TensorAlgebra}) = diffvars(Manifold(t))≠0
@pure isdyadic(::T) where T<:TensorAlgebra = isdyadic(T)
@pure isdual(::T) where T<:TensorAlgebra = isdual(T)
@pure istangent(::T) where T<:TensorAlgebra = istangent(T)

@pure isbasis(x) = false
@pure isdyadic(::Int) = false
@pure isdual(::Int) = false
@pure istangent(::Int) = false

@inline value_diff(m::T) where T<:TensorTerm = (v=value(m);istensor(v) ? v : m)

for T ∈ (:T,:(Type{T}))
    @eval @pure isbasis(::$T) where T<:TensorAlgebra = false
end
@pure UInt(m::T) where T<:TensorTerm = UInt(basis(m))

@pure hasconformal(V) = hasinf(V) && hasorigin(V)
@pure hasinf(::Int) = false
@pure hasinf(t::M) where M<:Manifold = hasinf(Manifold(t))
#@pure hasinf(m::T) where T<:TensorAlgebra = hasinf(Manifold(m))
@pure hasorigin(::Int) = false
@pure hasorigin(t::M) where M<:Manifold = hasorigin(Manifold(t))
#@pure hasorigin(m::T) where T<:TensorAlgebra = hasorigin(Manifold(m))

@pure hasorigin(V,B::UInt) = hasinf(V) ? (UInt(2)&B)==UInt(2) : isodd(B)

@pure hasinf(V,A::UInt,B::UInt) = hasconformal(V) && (isodd(A) || isodd(B))
@pure hasorigin(V,A::UInt,B::UInt) = hasconformal(V) && (hasorigin(V,A) || hasorigin(V,B))

@pure hasinf2(V,A::UInt,B::UInt) = hasconformal(V) && isodd(A) && isodd(B)
@pure hasorigin2(V,A::UInt,B::UInt) = hasconformal(V) && hasorigin(V,A) && hasorigin(V,B)

@pure function hasorigininf(V,A::UInt,B::UInt)
    hasconformal(V) && hasorigin(V,A) && isodd(B) && !hasorigin(V,B) && !isodd(A)
end
@pure function hasinforigin(V,A::UInt,B::UInt)
    hasconformal(V) && isodd(A) && hasorigin(V,B) && !isodd(B) && !hasorigin(V,A)
end

@pure hasinf2origin(V,A::UInt,B::UInt) = hasinf2(V,A,B) && hasorigin(V,A,B)
@pure hasorigin2inf(V,A::UInt,B::UInt) = hasorigin2(V,A,B) && hasinf(V,A,B)

@pure diffmask(::Int) = zero(UInt)
@pure function diffmask(V::M) where M<:Manifold
    d = diffvars(V)
    if isdyadic(V)
        v = ((one(UInt)<<d)-1)<<(mdims(V)-2d)
        w = ((one(UInt)<<d)-1)<<(mdims(V)-d)
        d<0 ? (typemax(UInt)-v,typemax(UInt)-w) : (v,w)
    else
        v = ((one(UInt)<<d)-1)<<(mdims(V)-d)
        d<0 ? typemax(UInt)-v : v
    end
end

@pure function symmetricsplit(V,a)
    sm,dm = symmetricmask(V,a),diffmask(V)
    isdyadic(V) ? (sm&dm[1],sm&dm[2]) : sm
end

@pure function symmetricmask(V,a)
    d = diffmask(V)
    a&(isdyadic(V) ? |(d...) : d)
end

@pure function symmetricmask(V,a,b)
    d = diffmask(V)
    D = isdyadic(V) ? |(d...) : d
    aD,bD = (a&D),(b&D)
    return a&~D, b&~D, aD|bD, aD&bD
end

@pure function diffcheck(V,A::UInt,B::UInt)
    d,db = diffvars(V),diffmask(V)
    v = isdyadic(V) ? db[1]|db[2] : db
    hi = hasinf2(V,A,B) && !hasorigin(V,A,B)
    ho = hasorigin2(V,A,B) && !hasinf(V,A,B)
    (hi || ho) || (d≠0 && count_ones(A&v)+count_ones(B&v)>diffmode(V))
end

## conversions

@pure function mixed(V,ibk::UInt)
    N,D,VC = mdims(V),diffvars(V),isdual(V)
    return if D≠0
        A,B = ibk&(UInt(1)<<(N-D)-1),ibk&diffmask(V)
        VC ? (A<<(N-D))|(B<<N) : A|(B<<(N-D))
    else
        VC ? ibk<<N : ibk
    end
end

@pure function combine(v,w,iak::UInt,ibk::UInt)
    (isdual(v) ≠ isdual(w)) && throw(error("$v and $w incompatible"))
    V,W = supermanifold(v),supermanifold(w)
    return if istangent(V)||istangent(W)
        gV,gW = (typeof(V)<:Int ? V : grade(V)),(typeof(W)<:Int ? W : grade(W))
        gras1,gras2 = iak&(UInt(1)<<gV-1),ibk&(UInt(1)<<gW-1)
        diffs = (iak&diffmask(W))|(ibk&diffmask(W))
        gras1|(gras2<<gV)|(diffs<<mdims(W)) # A|(B<<(N-D))
    else
        iak|(ibk<<mdims(V)) # ibk
    end
end

## functors

@pure loworder(N::Int) = N
@pure supermanifold(N::Int) = N

## adjoint parities

@pure parityreverse(G) = isodd(Int((G-1)*G/2))
@pure parityinvolute(G) = isodd(G)
@pure parityclifford(G) = parityreverse(G)⊻parityinvolute(G)
const parityconj = parityreverse

## reverse

@pure grade_basis(V::Int,B) = B&(one(UInt)<<V-1)
@pure grade_basis(V,B) = B&(one(UInt)<<grade(V)-1)
@pure grade(V,B::UInt) = count_ones(grade_basis(V,B))

## anti-reverse

@pure pseudograde(V::Int,B::UInt) = V-grade(V,B)
@pure pseudograde(V,B::UInt) = grade(V)-grade(V,B)

# comparison (special case for scalars)

Base.isless(a::T,b::S) where {T<:TensorTerm{V,0},S<:TensorTerm{W,0}} where {V,W} = isless(value(a),value(b))
Base.isless(a::T,b) where T<:TensorTerm{V,0} where V = isless(value(a),b)
Base.isless(a,b::T) where T<:TensorTerm{V,0} where V = isless(a,value(b))
Base.:<=(x::T,y::S) where {T<:TensorTerm{V,0},S<:TensorTerm{W,0}} where {V,W} = isless(x,y) | (x == y)
Base.:<=(x::T,y) where T<:TensorTerm{V,0} where V = isless(x,y) | (x == y)
Base.:<=(x,y::T) where T<:TensorTerm{V,0} where V = isless(x,y) | (x == y)

# operations

export ⊕, χ, gdims

"""
    χ(::TensorAlgebra)

Compute the Euler characteristic χ = ∑ₚ(-1)ᵖbₚ.
"""
χ(t::T) where T<:TensorAlgebra = (B=gdims(t);sum([B[t]*(-1)^t for t ∈ 1:length(B)]))
χ(t::T) where T<:TensorTerm = χ(Manifold(t),UInt(basis(t)),t)
@inline χ(V,b::UInt,t) = iszero(t) ? 0 : isodd(count_ones(symmetricmask(V,b,b)[1])) ? 1 : -1

function gdims(t::T) where T<:TensorTerm
    B,N = UInt(basis(t)),mdims(t)
    g = count_ones(symmetricmask(Manifold(t),B,B)[1])
    Variables{N+1,Int}([g==G ? abs(χ(t)) : 0 for G ∈ 0:N])
end
function gdims(t::T) where T<:TensorGraded{V,G} where {V,G}
    N = mdims(V)
    out = zeros(Variables{N+1,Int})
    ib = indexbasis(N,G)
    for k ∈ 1:length(ib)
        @inbounds t[k] ≠ 0 && (out[count_ones(symmetricmask(V,ib[k],ib[k])[1])+1] += 1)
    end
    return out
end

## set theory ∪,∩,⊆,⊇

@pure ∪(x::T) where T<:Manifold = x
@pure ∪(a::A,b::B,c::C...) where {A<:Manifold,B<:Manifold,C<:Manifold} = ∪(a∪b,c...)

@pure ∩(x::T) where T<:Manifold = x
@pure ∩(a::A,b::B,c::C...) where {A<:Manifold,B<:Manifold,C<:Manifold} = ∩(a∩b,c...)

## complement parity

@pure parityrighthodge(V::Int,B::Int,G,N=nothing) = isodd(V)⊻parityright(V,B,G,N)
@pure paritylefthodge(V::Int,B::Int,G,N) = (isodd(G) && iseven(N)) ⊻ parityrighthodge(V,B,G,N)
@pure parityright(V::Int,B::Int,G,N=nothing) = isodd(B+Int((G+1)*G/2))
@pure parityleft(V::Int,B::Int,G,N) = (isodd(G) && iseven(N)) ⊻ parityright(V,B,G,N)

for side ∈ (:left,:right)
    p = Symbol(:parity,side)
    pg = Symbol(p,:hodge)
    pn = Symbol(p,:null)
    pnp = Symbol(pn,:pre)
    @eval begin
        @pure $p(V::UInt,B::UInt,N::Int) = $p(0,sum(indices(B,N)),count_ones(B),N)
        @pure $pg(V::UInt,B::UInt,N::Int) = $pg(count_ones(V&B),sum(indices(B,N)),count_ones(B),N)
        @inline function $pn(V,B,v)
            if hasconformal(V) && count_ones(B&UInt(3)) == 1
                isodd(B) ? (2v) : (v/2)
            else v end
        end
        @inline function $pnp(V,B,v)
            if hasconformal(V) && count_ones(B&UInt(3)) == 1
                isodd(B) ? Expr(:call,:*,2,v) : Expr(:call,:/,v,2)
            else v end
        end
        #=@pure function $pg(V::Int,B::UInt,G=count_ones(B))
            ind = indices(B&(UInt(1)<<(mdims(V)-diffvars(V))-1),mdims(V))
            c = hasconformal(V) && (B&UInt(3) == UInt(2))
            $p(0,sum(ind),G,mdims(V)-diffvars(V))⊻c ? -1 : 1
        end=#
    end
end

@pure function complement(N::Int,B::UInt,D::Int=0,P::Int=0)::UInt
    UP,ND = UInt(1)<<(P==1 ? 0 : P)-1, N-D
    C = ((~B)&(UP⊻(UInt(1)<<ND-1)))|(B&(UP⊻((UInt(1)<<D-1)<<ND)))
    count_ones(C&UP)≠1 ? C⊻UP : C
end

# Hodge star ★

import AbstractTensors: complementrighthodge, complementright

## complement

export complementleft, complementright, ⋆, complementlefthodge, complementrighthodge

@doc """
    complementrighthodge(ω::TensorAlgebra)

Grassmann-Poincare-Hodge complement: ⋆ω = ω∗I
""" complementrighthodge

@doc """
    complementlefthodge(ω::TensorAlgebra)

Grassmann-Poincare left complement: ⋆'ω = I∗'ω
""" complementlefthodge

@doc """
    complementright(::TensorAlgebra)

Non-metric variant of Grassmann-Poincare-Hodge complement.
""" complementright

@doc """
    complementleft(::TensorAlgebra)

Non-metric variant Grassmann-Poincare left complement.
""" complementleft

# QR compatibility

@inline function LinearAlgebra.reflectorApply!(x::AbstractVector, τ::TensorAlgebra, A::StridedMatrix)
    @assert !LinearAlgebra.has_offset_axes(x)
    m, n = size(A)
    if length(x) != m
        throw(DimensionMismatch("reflector has length $(length(x)), which must match the dimension of matrix A, $m"))
    end
    @inbounds begin
        for j = 1:n
            #dot
            vAj = A[1,j]
            for i = 2:m
                vAj += Base.conj(x[i])*A[i,j]
            end

            vAj = conj(τ)*vAj

            #ger
            A[1, j] -= vAj
            for i = 2:m
                A[i,j] -= x[i]*vAj
            end
        end
    end
    return A
end
