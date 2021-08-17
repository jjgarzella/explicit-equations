module systemdensity

using HomotopyContinuation
using Oscar
using Nemo
using Singular

function count_nonzero_coeff(S,vars)
    t = length(vars)
    T = t+binomial(t,t-1)+binomial(t+1,t-1)
    U = Set()

    for f in S
        countf = HomotopyContinuation.monomials(f)
        for mon in countf
            push!(U,mon)
        end
    end
    V = length(U)
    println(V/T)
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
"""



end #module systemdensity
