
#   This file is part of Leibniz.jl. It is licensed under the AGPL license
#   Leibniz Copyright (C) 2019 Michael Reed

# alpha-numeric digits
const vio = ('∞','∅')
const digs = "1234567890"
const low_case,upp_case = "abcdefghijklmnopqrstuvwxyz","ABCDEFGHIJKLMNOPQRSTUVWXYZ"
const low_greek,upp_greek = "αβγδϵζηθικλμνξοπρστυφχψω","ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΡΣΤΥΦΨΩ"
const alphanumv = digs*low_case*upp_case #*low_greek*upp_greek
const alphanumw = digs*upp_case*low_case #*upp_greek*low_greek

# subscript index
const subs = Dict{Int,Char}(
   -1 => vio[1],
    0 => vio[2],
    1 => '₁',
    2 => '₂',
    3 => '₃',
    4 => '₄',
    5 => '₅',
    6 => '₆',
    7 => '₇',
    8 => '₈',
    9 => '₉',
    10 => '₀',
    [j=>alphanumv[j] for j ∈ 11:36]...
)

# superscript index
const sups = Dict{Int,Char}(
   -1 => vio[1],
    0 => vio[2],
    1 => '¹',
    2 => '²',
    3 => '³',
    4 => '⁴',
    5 => '⁵',
    6 => '⁶',
    7 => '⁷',
    8 => '⁸',
    9 => '⁹',
    10 => '⁰',
    [j=>alphanumw[j] for j ∈ 11:36]...
)

# vector and co-vector prefix
const pre = ("v","w","∂","ϵ")
const PRE = ("X","x","Y","y")

# vector space and dual-space symbols
const vsn = (:V,:VV,:W)
const VSN = (:Χ,:ΧΧ,:Υ) # \Chi,\Upsilon

# converts indices into BitArray of length N
@inline function indexbits(N::I,indices::T) where {I<:Integer,T<:SVTI}
    out = falses(N)
    for k ∈ indices
        out[k] = true
    end
    return out
end

# index sets
index_limit = 20
const digits_fast_cache = Vector{Values}[]
const digits_fast_extra = Dict{UInt,Values}[]
@pure digits_fast_calc(b,N) = Values{N+1,Int}(digits(b,base=2,pad=N+1))
@pure function digits_fast(b,N)
    if N>index_limit
        n = N-index_limit
        for k ∈ length(digits_fast_extra)+1:n
            push!(digits_fast_extra,Dict{UInt,Values{k+1,Int}}())
        end
        !haskey(digits_fast_extra[n],b) && push!(digits_fast_extra[n],b=>digits_fast_calc(b,N))
        @inbounds digits_fast_extra[n][b]
    elseif N==0
        Values{1,Int}(0)
    else
        for k ∈ length(digits_fast_cache)+1:min(N,index_limit)
            push!(digits_fast_cache,[digits_fast_calc(d,k) for d ∈ 0:1<<(k+1)-1])
            GC.gc()
        end
        @inbounds digits_fast_cache[N][b+1]
    end
end

const indices_cache = Dict{UInt,Vector{Int}}()
indices(b::UInt) = findall(digits(b,base=2).==1)
function indices_calc(b::UInt,N::Int)
    d = digits_fast(b,N)
    l = length(d)
    a = Int[]
    for i ∈ 1:l
        d[i] == 1 && push!(a,i)
    end
    return a
end
function indices(b::UInt,N::Int)
    !haskey(indices_cache,b) && push!(indices_cache,b=>indices_calc(b,N))
    return @inbounds indices_cache[b]
end

shift_indices(V::T,b::UInt) where T<:Manifold = shift_indices(supermanifold(V),b)
function shift_indices!(s::T,set::Vector{Int}) where T<:Manifold
    M = supermanifold(s)
    if !isempty(set)
        k = 1
        hasinf(M) && set[1] == 1 && (set[1] = -1; k += 1)
        shift = hasinf(M) + hasorigin(M)
        hasorigin(M) && length(set)>=k && set[k]==shift && (set[k]=0;k+=1)
        shift > 0 && (set[k:end] .-= shift)
    end
    return set
end

shift_indices(V::Int,b::UInt) = indices(b,V)
shift_indices!(s::Int,set::Vector{Int}) = set

# printing of indices
@inline function printindex(i,l::Bool=false,e::String=pre[1],pre=pre)
    t = i>36; j = t ? i-26 : i
    (l&&(0<j≤10)) ? j : ((e∉pre[[1,3]])⊻t ? sups[j] : subs[j])
 end
@inline printindices(io::IO,b::UInt,l::Bool=false,e::String=pre[1],pre::NTuple{4,String}=pre) = printindices(io,indices(b),l,e,pre)
@inline printindices(io::IO,b::VTI,l::Bool=false,e::String=pre[1],pre::NTuple{4,String}=pre) = print(io,e,[printindex(i,l,e,pre) for i ∈ b]...)
@inline printindices(io::IO,a::VTI,b::VTI,l::Bool=false,e::String=pre[1],f::String=pre[2]) = printindices(io,a,b,Int[],Int[],l,e,f)
@inline function printindices(io::IO,a::VTI,b::VTI,c::VTI,d::VTI,l::Bool=false,e::String=pre[1],f::String=pre[2],g::String=pre[3],h::String=pre[4])
    A,B,C,D = isempty(a),!isempty(b),!isempty(c),!isempty(d)
    PRE = (e,f,g,h)
    C && printindices(io,c,l,g,PRE)
    D && printindices(io,d,l,h,PRE)
    !((B || C || D) && A) && printindices(io,a,l,e,PRE)
    B && printindices(io,b,l,f,PRE)
end

printindices(io::IO,V::Int,e::UInt,label::Bool=false) = printlabel(io,V,e,label,pre...)
@inline function printlabel(io::IO,V::Int,e::UInt,label::Bool,vec,cov,duo,dif)
    printindices(io,shift_indices(V,e),label,vec)
    return io
end

@inline function printlabel(io::IO,V::T,e::UInt,label::Bool,vec,cov,duo,dif) where T<:Manifold
    M = supermanifold(V)
    N,D,C,db = mdims(M),diffvars(M),dyadmode(V),diffmask(V)
    if C < 0
        es = e & (~(db[1]|db[2]))
        n = Int((N-2D)/2)
        eps = shift_indices(V,e & db[1]).-(N-2D-hasinf(M)-hasorigin(M))
        par = shift_indices(V,e & db[2]).-(N-D-hasinf(M)-hasorigin(M))
        printindices(io,shift_indices(V,es & UInt(2^n-1)),shift_indices(V,es>>n),eps,par,label,vec,cov,duo,dif)
    else
        es = e & (~db)
        eps = shift_indices(V,e & db).-(N-D-hasinf(M)-hasorigin(M))
        if !isempty(eps)
            printindices(io,shift_indices(V,es),Int[],C>0 ? Int[] : eps,C>0 ? eps : Int[],label,C>0 ? cov : vec,cov,C>0 ? dif : duo,dif)
        else
            printindices(io,shift_indices(V,es),label,C>0 ? string(cov) : vec)
        end
    end
    return io
end

@inline printlabel(V::T,e::UInt,label::Bool,vec,cov,duo,dif) where T<:Manifold = printlabel(IOBuffer(),V,e,label,vec,cov,duo,dif) |> take! |> String

showparens(T) = !|(broadcast(<:,T,parnot)...) && |(broadcast(<:,T,parval)...)

function showstar(io::IO,v)
    if !(isa(v,Integer) && !isa(v,Bool) || isa(v,AbstractFloat) && isfinite(v))
        print(io, "*")
    end
end

function showvalue(io::IO,V,B::UInt,i::T) where T
    if showparens(T)
        print(io,"(",i,")")
    else
        show(io,i)
        showstar(io,i)
    end
    printindices(io,V,B)
end

function indexstring(V::M,D) where M<:Manifold
    io = IOBuffer()
    printlabel(io,V,D,true,PRE...)
    String(take!(io))
end

indexsymbol(V::M,D) where M<:Manifold = Symbol(indexstring(V,D))

indexsplit(B,N) = [UInt(1)<<(k-1) for k ∈ indices(B,N)]

indexparity!(ind::Values{N,Int}) where N = indexparity!(Variables(ind))
function indexparity!(ind::Variables{N,Int}) where N
    k = 1
    t = false
    while k < length(ind)
        if ind[k] > ind[k+1]
            ind[Values(k,k+1)] = ind[Values(k+1,k)]
            t = !t
            k ≠ 1 && (k -= 1)
        else
            k += 1
        end
    end
    return t, ind
end
function indexparity!(ind::Vector{Int},s)
    k = 1
    t = false
    while k < length(ind)
        if ind[k] == ind[k+1]
            ind[k] == 1 && hasinf(s) && (return t, ind, true)
            isone(s[ind[k]]) && (t = !t)
            deleteat!(ind,[k,k+1])
        elseif ind[k] > ind[k+1]
            ind[[k,k+1]] = ind[[k+1,k]]
            t = !t
            k ≠ 1 && (k -= 1)
        else
            k += 1
        end
    end
    return t, ind, false
end
