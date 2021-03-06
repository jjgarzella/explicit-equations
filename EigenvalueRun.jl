module EigenvalueRun

using DynamicPolynomials
using LinearAlgebra
using GenericSchur
using GenericSVD
using Statistics
using SmithNormalForm  # to add: pkg> add https://github.com/wildart/SmithNormalForm.jl.git

import TaylorSeries

include("JuliaEigenvalueSolver/src/EigenvalueSolver.jl")

@polyvar c[1:20]
@polyvar d[1:6]

x = vcat(c,d)

function run()
    S = Polynomial{true, ComplexF64}[]
    for n = 1:67
        S = vcat(S,convert(Vector{Polynomial{true,ComplexF64}},eval(Meta.parse(readline("separatedBFrels/RelsEq$n")))))
    end
    S
    solv_dense = @timed EigenvalueSolver.solve_OD_dense(S, x; maxdeg = 25, complex = true, verbose=false, trySemiregular=false)
    return solv_dense
    return solv_dense[1]
end

"""
@polyvar x[1:3]
H = [x[1]^2+x[3]^2, x[1]^2+(1+sqrt(complex(-15)))*x[2]-8*x[3],x[1]^2+x[2]^2+x[3]^3,-x[1]^2+x[3]^4-8]
solv_dense = @timed EigenvalueSolver.solve_OD_dense(H, x; maxdeg = 15, complex = false, verbose=false, trySemiregular=false)

W = [x[1]^2-1,x[1]-1]
solv_dense = @timed EigenvalueSolver.solve_OD_dense(W, x; maxdeg = 25, complex = false, verbose=false, trySemiregular=false)

G = [x[1]^2-1,x[1]-1,x[1]*x[2]+1]
solv_dense = @timed EigenvalueSolver.solve_OD_dense(G, x; maxdeg = 25, complex = false, verbose=false, trySemiregular=false)


# Note to Kevin
# include("JuliaEigenvalueSolver/src/EigenvalueSolver.jl")
# DynamicPolynomials is similar to HomotopyContinuation in
# polynomial input.
# Test @polyvar x[1:3]
# H = [x[1]^2+x[3]^2, x[1]^2+(1+sqrt(complex(-15)))*x[2]-8*x[3],x[1]^2+x[2]^2+x[3]^3,-x[1]^2+x[3]^4-8]

"""




end
