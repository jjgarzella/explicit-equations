module Homotopy3Run

using HomotopyContinuation
using DynamicPolynomials

@polyvar c[1:20]
@polyvar d[1:6]

function setup()
    S = Polynomial{true, ComplexF32}[]
    for n = 1:67
        S = vcat(S,convert(Vector{Polynomial{true,ComplexF64}},eval(Meta.parse(readline("separatedBFrels/RelsEq$n")))))
    end
    S
end

function run()
    S = Homotopy3Run.setup()
    #H = [d[1]^2,d[2]+d[1]-d[3]+5,d[2]-3,d[1]^2*d[1]-d[2]*d[3]]
    @time HomotopyContinuation.solve(S,threading = true)
end


end #module
