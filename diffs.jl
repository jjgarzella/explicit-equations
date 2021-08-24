module diffs

using HomotopyContinuation
using Oscar
using Nemo
using Singular
using StatsBase

include("VMLS.jl/src/VMLS.jl")

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

# the Borisov trick is in chooserand

function chooserand(S,point,i,allvars)
    vecs = []
    j = 1
    A = sample(1:length(S),i,replace=false)
        for a in A
            push!(vecs,Array{Int64}(diffs.differentials(S[a],point,allvars)))
        end
    B=vecs
    while length(VMLS.gram_schmidt(B)) !== length(B) && j <= binomial(length(S),i)
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

    println(A)
    return T
end
# The Array{Int64} is needed to remove an error that arises with
# the gram_schmidt algorithm

# a more "meticulous" application of Borisov's trick

function choosemetic(S,point,i,j,allvars)
    X = 1
    vecs = []
    vecs2 = []
    while X == 1
        B = sample(1:length(S),j,replace=false)
        for b in B
            push!(vecs2,S[b])
        end
        j = 1
        A = sample(1:length(vecs2),i,replace=false)
        for a in A
            push!(vecs,Array{Int64}(diffs.differentials(vecs2[a],point,allvars)))
        end
        C = vecs

        if length(VMLS.gram_schmidt(C)) == length(C)
            X = 2
            T = []
            for a in A
                push!(T,S[a])
            end
            println(A)
            return T
        elseif length(VMLS.gram_schmidt(C)) !== length(C)
            D = []
            while length(VMLS.gram_schmidt(C)) !== length(C) && j <= binomial(length(B),i)
                D = sample(1:length(B),i,replace=false)
                j += 1
                vecs = []
                for d in D
                    push!(vecs,Array{Int64}(diffs.differentials(S[d],point,allvars)))
                end
                C = vecs
                if j == length(vecs2)
                    println("All possible choices of $i equations from the system yielded linearly
dependent differentials at the point.")
                end
            end
            println(D)
            println(C)
            X = 2
            T = []
            for d in D
                push!(T,S[d])
            end
            return T
        end
    end
end

"""
Example 3:

@polyvar x[1:3]
allvars = x
S = [x[1]^2+x[2]^2, x[1]^2+x[3]^2, x[2]^2+x[1]x[3], x[2]^2+x[1]x[2]x[3]]
point = [1,4,5]
diffs.chooserand(S,point,2,allvars)

Examples 4:

@polyvar x[1:2]
allvars = x
S = [x[1]+x[2]^2, x[1]^2+2*x[2], x[1]^2, x[2]+2*x[1]^2,x[1]-x[2]+x[2]^2]
point = [1/2,1]
diffs.chooserand(S,point,3,allvars)
diffs.choosemetic(S,point,2,5,allvars)
"""

"""
To run on Borisov System, we will need to do the following.

set up the equations
input the point we are extending
then write in diffs.choosemetic(S,point,16,32,allvars)
"""


end # module diffs
