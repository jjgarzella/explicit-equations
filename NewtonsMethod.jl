module NewtonsMethod

using HomotopyContinuation
using Oscar
using Nemo
using Singular
using GroebnerBasis
using DynamicPolynomials

@polyvar c[1:20]
@polyvar d[1:6]

vars = vcat(c,d)

function inputBF()
  S = DynamicPolynomials.Polynomial{true, BigInt}[]
  for n = 1:67
      S = vcat(S,convert(Vector{DynamicPolynomials.Polynomial{true,BigInt}},eval(Meta.parse(readline("separatedBFrels/RelsEqMod$n")))))
  end
  S
end

function run(S,p)
  @polyvar c[1:20]
  @polyvar d[1:6]
  vars = vcat(c,d)

  newtons_method_p(S,p,vars)
end

#please run this one for the computations

function setupandrun(p)
    #S = DynamicPolynomials.Polynomial{true, ComplexF64}[]
    #for n = 1:67
    #    S = vcat(S,convert(Vector{DynamicPolynomials.Polynomial{true,ComplexF64}},eval(Meta.parse(readline("separatedBFrels/RelsEqMod$n")))))
    #end
    #S

    S = DynamicPolynomials.Polynomial{true, BigInt}[]
    for n = 1:67
        S = vcat(S,convert(Vector{DynamicPolynomials.Polynomial{true,BigInt}},eval(Meta.parse(readline("separatedBFrels/RelsEqMod$n")))))
    end
    #Int128 gives inexact error
    #S = Nemo.gfp_mpoly[]
    #for n = 1:67
    #    S = vcat(S,convert(Vector{Nemo.gfp_mpoly},eval(Meta.parse(readline("separatedBFrels/RelsEqMod$n")))))
    #end
    #S
    println(S)

    newtons_method_p(S,p,vars)
end


"""
@polyvar x[1:2]
var = x
p = 3
f = DynamicPolynomials.Polynomial{true,ComplexF64}[223*x[1]^3-2x[1]*x[2],x[1]^2*x[2]-2x[2]^2+x[1]]
NewtonsMethod.newtons_method_p(f,p,var)

@polyvar x[1:2]
var = x
p = 3
f = [x[1]-x[2],x[2]]
NewtonsMethod.newtons_method_p(f,p,var)


@polyvar x[1:2]
@polyvar y[1:2]
var = vcat(x,y)
p = 3
f = [x[1]-x[2],x[2]-y[1]+x[1]+y[2],y[2]+y[1]-x[1]+x[2],x[1]^2+y[1]-y[2]+x[2],x[2]+y[2]-y[2]]
NewtonsMethod.newtons_method_p(f,p,var)

@polyvar x[1:3]
var = x
p = 5
f = [x[1]-x[2],x[1]^2-x[3],x[2],x[2]^2-3*x[3]]
NewtonsMethod.newtons_method_p(f,p,var)

@polyvar x[1:3]
var = x
p = 7
f = [x[1]*x[3] - x[2]^2, x[2] - x[3]^2, x[1] - x[2]*x[3], x[1] + x[2] + x[3] + 1]
NewtonsMethod.newtons_method_p(f,p,var)

"""


function newtons_method_p(f,p,var)
  sol_modp =   NewtonsMethod.convertarraygfptoint(NewtonsMethod.solve_modwithrat(f,p,var)[1])
  #NewtonsMethod.solve_modwithrat(f,p,var)[1]
  sol_mod_pn = sol_modp
  println(sol_modp)
  for n in 2:200
    sol_mod_pn = NewtonsMethod.extend_mod(sol_mod_pn,f,var,n,p)
    #println("$n: ")
    #println(NewtonsMethod.extend_mod(sol_mod_pn,f,var,n,p))
    println(sol_mod_pn)
  end
  final_sol = sol_mod_pn
  exact = find_exact_sol(final_sol,[sqrt(complex(-15)),p^200])
  println(exact)
  if double_check(f,final_sol,var) == true
    println("The final solution solves the system mod $p raised 200. The final solution
and its exact form are respectively in the array:")
  return [final_sol,exact]
  else
    println("The final solution obtained does not solve the system in mod $p raised 5.")
    return [final_sol,exact]
  end
end

function solvegroebnerp(I,p,var)
  G = gens(I)
  sol = zeros(Int,length(var))
  for i = 1:min(length(G),length(var)) # maybe just up to length(var)
    Gvars = vars(G[i])
    indices = Int[]

    for j = 1:length(Gvars)
      push!(indices,findfirst(isequal(Gvars[j]),var))
    end

    for k = 1:p
      if i <= length(indices)
        sol[indices[i]] = k
        if Singular.evaluate(G[i],sol)== 0
          sol[indices[i]] = mod(k,p)
        elseif k == p && Singular.evaluate(G[i],sol) !== 0
          sol[indices[i]] = 0
          println("Error: Could not solve equation $i within mod p.")
        else
          nothing
        end
      elseif i > length(indices) && k == 1
        println("Error: Indexing failure in function solvegroebnerp.")
      else
        nothing
      end
    end
  end
sol
end#solvegroeber

function makesystem(f,p,var)
  F = Fp(p)
  variablesnew = string.(var)
  R, (indeter) = Singular.PolynomialRing(F,variablesnew)
  #println(Singular.PolynomialRing(F,variablesnew))
  g = repeat([indeter[1]],length(f))
  for i in 1:length(f)
    g[i] = f[i]((var...,) => (indeter...,))
  end
  return g
end

function solve_mod(f,p,var)
  F = Fp(p)
  variablesnew = string.(var)
  R, (indeter) = Singular.PolynomialRing(F,variablesnew)
  #println(Singular.PolynomialRing(F,variablesnew))
  R = makering(f,p,var)
  g = repeat([indeter[1]],length(f))
  for i in 1:length(f)
    g[i] = f[i]((var...,) => (indeter...,))
  end

  I = Ideal(R,g)

  I_groebner = GroebnerBasis.f4(I)
  # Need to solve mod p
  print(I_groebner)
  print("
  WARNING: Make sure to check that the Groebner Basis is solvable via function method of solvegroebnerp.
  In particular, check that first entry consists of only a single variable and the nth entry has
  only n distinct variables.
  Solution may not be accurate if not.
")
  solvegroebnerp(I,p,indeter)
end

"""
Example for solve_mod type in @polyvar x[1:3], f = [x[1]^2+1,x[2]+1,x[1]]
p=3, and var = x.
the program should fail since this system cannot be solved

for another test, try f=[x[1]^3-x[1]^2+2, x[2]-2x[1]^2+1, x[3]-3x[1]+5]
the Groebner basis is (x[3] - 1, x[1]^2 + x[2] + 1, x[1]*x[2] + x[1] - x[2], x[2]^2 - x[1] + 1)
and the program fails to solve this

for a currently successful test, try f = [x[1]^3-2x[1]*x[2],x[1]^2*x[2]-2x[2]^2+x[1]]
which gives Groebner basis (x[1]^2, x[2]^2 + 2*x[1], x[1]*x[2])
and solution [0,0].
"""

function solve_modwithrat(f,p,var)
  F = GF(p)
  variablesnew = string.(var)
  R, (indeter) = Oscar.PolynomialRing(F,variablesnew)
  println(indeter)
  g = [indeter[1]]

  for n in 1:length(f)-1
    push!(g,indeter[1])
  end
  #println(g)
  #g = repeat([indeter[1]],length(f))
  #println(g)
  #println(f)
  for i in 1:length(f)
    g[i] = f[i]((var...,) => (indeter...,))
    # g[i] = HomotopyContinuation.evaluate(f[i], vars => indeter)
    println(g[i])
  end

  I = Oscar.ideal(g)
  #println(I)
  if !isequal(points(I),[])
    sol = points(I)
  elseif !isequal(proj_points(I),[])
    sol = proj_points(I)
  else
    sol = [zeros(Int64,length(var))]
    println("There was an error with rational_points. Potentially, no solutions were found.")
  end
  sol
end

function convertgfptofloat(gfp)
  parse(Float64,string(gfp))
end

function convertgfptoint(gfp)
  parse(Int64,string(gfp))
end

function convertarraygfptoint(array)
  newarray = zeros(Int,length(array))
  for i = 1:length(array)
    newarray[i] =convertgfptoint(array[i])
  end
  newarray
end

function pseudo_inv_finite_field(M,p)
  # input M for matrix
  # formula from wikipedia: A^+ = (A^*A)^{-1}A^*
  # or A^*(A^*A)^{-1} depending on whether A has LI col/row res.
  if rank(M) == size(M)[1]
    pseudoinv = mod.(mod.(transpose(M),p)*mod.(inv(M*transpose(M)),p),p)
  elseif rank(M) == size(M)[2]
    pseudoinv = mod.(mod.(inv(transpose(M) * M),p) * mod.(transpose(M),p),p)
  else
    println("Cannot compute pseudoinverse of matrix input. Check rank.")
  end

  # only if our matrices are full rank is this well-defined
  # example: F[1 2; 2 3] gives 2x2 matrix with entries 1,2 ; 2,0.
  # Current issue is that M needs to be defined via above way before inputting.
  # this goes not compute M over modulo p. One must construct
  # M similarly to F[1 2; 2 4] by setting F = GF(p) for p prime.
  # Needs HomotopyContinuation to run.
  pseudoinv = round.(pseudoinv)

  pseudoinv = convert(Array{Int},pseudoinv)
  return pseudoinv
end


function jacobian_sol(f,sol,p)
  # makes Jacobian matrix of system f of polynomials and evalutes at sol
  # example: f =[x[1]^2+1,x[2]*x[1]]
  F = System(f)
  dummyjacob = jacobian(F,sol)
  for i in 1:nrows(dummyjacob)
    for j in 1:ncols(dummyjacob)
      dummyjacob[i,j] = dummyjacob[i,j] % p
    end
  end
  return dummyjacob
end

function evaluatesystem(f,sol,var)
  check(g) = double_check_single(g,sol,var)
  # broadcasts the function to the array
  output = check.(f)
  return output
end

function extend_mod(sol,f,var,m,p)
  # a_{n+1} = a_n - Df(a_n)^+ * f(a_n)
  # println("Extending mod solution mod $p raised $m.")
  for n = 1:m
    extendsol = sol - NewtonsMethod.pseudo_inv_finite_field(NewtonsMethod.jacobian_sol(f,sol,p),p)*NewtonsMethod.evaluatesystem(f,sol,var)
    #println("Extension $n:", extendsol)
    sol = extendsol
  end
  for n = 1:length(sol)
    sol[n] = mod(sol[n],p^m)
  end
  return sol

  # test to check this works. Run in repl the following.
  # @polyvar x[1:3], f = [x[1]+1, x[1]x[2]+2, x[1]+x[2]+x[3]],
  # sol = [1, 3, 5], NewtonsMethod.extend_mod(sol,f,x,10,3).
  # observe sol does not satisfy f. however, the output should be [-1, 2,-1] mod 3^(10)
  # which is then [59048.0,2.0,59048.0]
  # which does satisfy f mod 3^(10)
end

function powersof(a,n)
  # puts powers of a into an array (up to n)
  A = [1]
  for i = 1:n
    push!(A,a^i)
  end
  A
end

function find_exact_sol(sol,addrelation)
  # addrelation; W = [sqrt(complex(-15)),ppower]
  # lindep (in Nemo)
  # p is the prime or prime power
  CC = ComplexField(128) # 64-bits of accuracy
  V = []
  for a in sol
    U = [CC(1)]
    push!(U,CC(real(a),imag(a)))
    for i in addrelation
    push!(U,CC(real(i),imag(i)))
    end
    push!(V,Nemo.lindep(U,20))
  end

  EV=[]

  for b in  V
    esol =
    complex(((-NewtonsMethod.convertgfptofloat(b[1]))/-NewtonsMethod.convertgfptofloat(b[2])) -
    (-NewtonsMethod.convertgfptofloat(b[3])*sqrt(complex(-15)))/(-NewtonsMethod.convertgfptofloat(b[2])))

    push!(EV,esol)
  end

  EV
end

function double_check_single(g,sol,var)
  #println(g)
  # for n in 1:length(sol)
  #   g=subs(g, (vs...,)=>(sol...,))
  # end

  g( (var...,) => (sol...,))
end

function double_check(f,sol,var)
  # output = []
  # for g in f
  #   double_check_single(g,sol)
  #   push!(output, double_check_single(g,sol))
  # end
  # output

  # defines a function
  check(g) = double_check_single(g,sol,var)
  # broadcasts the function to the array
  output = check.(f)
  # println(output)
  all(output .== zeros(length(output)))
end


#rationalpoints functions


# evaluate function for an ideal
function idealevaluate(I::Oscar.MPolyIdeal, v)
    R = base_ring(I)
    return ideal([R(Oscar.evaluate(f,v)) for f in gens(I)])
end

# convert to univariate
function to_univar(a, var::Int, Kx)
    r = zero(Kx)
    for (c, ev) in zip(Oscar.coefficients(a), exponent_vectors(a))
        setcoeff!(r, ev[var], c)
    end
    return r
end

# enumeration of points; starts with a partial solution `part`
# unsafe: points are pushed into `ans`
# would be nice to have multithread
function enum!(I, part, els, ans)
    R = base_ring(I)
    char = characteristic(R)
    F = base_ring(R)
    n = nvars(R)
    i = length(part)
    if i == n
        if I == ideal([R(0)]) push!(ans, part) end
        return
    end
    rest = gens(R)[i+2:end]
    if i < n-1
        elim = gens(Oscar.eliminate(I, rest))[1]
    else
        # eliminate an empty list doesn't work right now
        I = ideal(groebner_basis(I))
        elim = gens(I)[1]
    end
    if elim == 0
        values = els
    else
        values = roots(to_univar(elim, i+1, F["x"][1]))
    end
    for v in values
        vec = [R(u) for u in vcat(part, [v], rest)]
        enum!(idealevaluate(I, vec), vcat(part, [v]), els, ans)
    end
end

# enumeration of points for an affine scheme
function points(I::Oscar.MPolyIdeal)
    els = collect(base_ring(base_ring(I)))
    ans = []
    enum!(I, [], els, ans)
    ans
end

# enumeration of points for a projective scheme
function proj_points(I::Oscar.MPolyIdeal)
    R = base_ring(I)
    F = base_ring(R)
    n = nvars(R)
    els = collect(F)
    ans = []
    for i in 1:n
        rest = gens(R)[i+1:end]
        part = vcat(repeat([F(0)], i-1), [F(1)])
        vec = [R(u) for u in vcat(part, rest)]
        enum!(idealevaluate(I, vec), part, els, ans)
    end
    ans
end

# example taken from the documentation of Magma
# function example(char=7823)
#     R,(x,y,z,w) = Oscar.PolynomialRing(GF(char),["x","y","z","w"])
#     I = ideal([4*x*z+2*x*w+y^2+4*y*w+7821*z^2+7820*w^2,4*x^2+4*x*y+7821*x*w+7822*y^2+7821*y*w+7821*z^2+7819*z*w+7820*w^2])
#     @time length(proj_points(I))
# end
# over an extension field

# function exampleGF()
#     R,(x,y) = PolynomialRing(FiniteField(19,3,"a")[1],["x","y"])
#     I = ideal([y^27-x^8*(1-x)])
#     @time length(points(I))
# end






end#module
