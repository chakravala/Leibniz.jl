
#   This file is part of Leibniz.jl. It is licensed under the AGPL license
#   Leibniz Copyright (C) 2019 Michael Reed

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
            else
                v
            end
        end
        @inline function $pnp(V,B,v)
            if hasconformal(V) && count_ones(B&UInt(3)) == 1
                isodd(B) ? Expr(:call,:*,2,v) : Expr(:call,:/,v,2)
            else
                v
            end
        end
        @pure function $pg(V::Int,B::UInt,G=count_ones(B))
            ind = indices(B&(UInt(1)<<(mdims(V)-diffvars(V))-1),mdims(V))
            c = hasconformal(V) && (B&UInt(3) == UInt(2))
            $p(0,sum(ind),G,mdims(V)-diffvars(V))⊻c ? -1 : 1
        end
    end
end

@pure function complement(N::Int,B::UInt,D::Int=0,P::Int=0)::UInt
    UP,ND = UInt(1)<<(P==1 ? 0 : P)-1, N-D
    C = ((~B)&(UP⊻(UInt(1)<<ND-1)))|(B&(UP⊻((UInt(1)<<D-1)<<ND)))
    count_ones(C&UP)≠1 ? C⊻UP : C
end

# Hodge star ★

const complementrighthodge = ⋆
const complementright = !

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
