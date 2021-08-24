module SystemDensityRun

using HomotopyContinuation
using Oscar
using Nemo
using Singular

include("JuliaEigenvalueSolver/src/EigenvalueSolver.jl")

function count_nonzero_coeff(S,vars)
    t = length(vars)
    T = 1+t+binomial(t+1,2)+binomial(t+2,3)
    U = Set()

    for f in S
        countf = HomotopyContinuation.monomials(f)
        for mon in countf
            push!(U,mon)
        end
    end
    V = length(U)
    println(V/T)
    println(V)
    println(T)
    return U
end

"""
Example 1:

@polyvar x[1:3]
S = [x[1]^2]
vars = x
systemdensity.count_nonzero_coeff(S,vars)

Example 2:

@polyvar x[1:3]
vars = x
S = [x[1]^2+x[2]^2, x[1]^2+x[3]^2, x[2]^2+x[1]x[3], x[2]^2+x[1]x[2]x[3],x[2]x[3]+x[2]]
systemdensity.count_nonzero_coeff(S,vars)

EigenvalueSolver.getRandomSystem_dense(x,[2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3];complex=true)
S = EigenvalueSolver.getRandomSystem_dense(x,[3,2,3,3,3,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3];complex=true)
SystemDensityRun.count_nonzero_coeff(S,x)
"""



end #module systemdensity
