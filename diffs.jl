module diffs

using HomotopyContinuation
using Oscar
using Nemo
using Singular
using StatsBase

include("VMLS/src/VMLS.jl")

# script to compute the differential of a polynomial at a point
# script to checck linear independence of a set of vector

function diff(f,vars,point,allvars)
    # all vars include all variables initialized using @var and @polyvar
    h = differentiate(f, vars)
    subs(h,allvars=>point)
end

"""
Example 1:

@polyvar x[1:3]
allvars = x
f = x[1]^2
vars = x[1]
point = [1,2,3]
diffs.diff(f,vars,point,allvars)
"""

function differentials(f,point,allvars)
    A = []
    for v in allvars
        push!(A,diff(f,v,point,allvars))
    end
    return A
end

"""
Example 2:

@polyvar x[1:3]
allvars = x
f = x[1]^2
vars = x[1]
point = [1,2,3]
diffs.differentials(f,point,allvars)
"""

function choose(S,point,i,allvars)
    vecs = []
    j = 1
    A = sample(1:length(S),i,replace=false)
        for a in A
            push!(vecs,Array{Int64}(diffs.differentials(S[a],point,allvars)))
        end
    B=vecs
    while length(VMLS.gram_schmidt(B)) !== length(B) & j <= binomial(length(S),i)
        A = sample(1:length(S),i,replace=false)
        j += 1
        vecs = []
            for a in A
                push!(vecs,Array{Int64}(diffs.differentials(S[a],point,allvars)))
            end
        B=vecs
        if j == length(S)
            println("All possible choices of $i equations from the system yielded linearly
            dependent differentials at the point.")
        end
    end

    T = []

    for a in A
        push!(T,S[a])
    end

    return T
end
# The Array{Int64} is needed to remove an error that arisesw with
# the gram_schmidt algorith.

"""
Example 3:

@polyvar x[1:3]
allvars = x
S = [x[1]^2+x[2]^2, x[1]^2+x[3]^2, x[2]^2+x[1]x[3], x[2]^2+x[1]x[2]x[3]]
point = [1,4,5]
diffs.choose(S,point,2,allvars)

Examples 4:

@polyvar x[1:2]
allvars = x
S = [x[1]+x[2]^2, x[1]^2+2*x[2], x[1]^2]
point = [1/2,1]
diffs.choose(S,point,3,allvars)
"""



end # module diffs
