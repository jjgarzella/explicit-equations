module test

using HomotopyContinuation
using Oscar
using Nemo
using Singular
using GroebnerBasis

include("JuliaEigenvalueSolver/src/EigenvalueSolver.jl")


function testmemoryii(k,m)
    for l = 1:k
        F = randomdenseii(m,l)
        println(F)
        println(@time HCsolve(F))
        println("The system of size $l has been successful.")
    end
end


function randomdenseii(n,m)
  @polyvar x[1:n]
  ds = []
  for i = 1:m
    push!(ds,i+1)
  end
  @time EigenvalueSolver.getRandomSystem_dense(x, ds; complex = true)
end

function testmemory(k)
    for l = 1:k
        F = randomdense(l,l)
        println(@time HCsolve(F))
        println("The system of size $l has been successful.")
    end
end


function randomdense(n,m)
  @polyvar x[1:n]
  ds = []
  for i = 1:m
    push!(ds,i)
  end
  @time EigenvalueSolver.getRandomSystem_dense(x, ds; complex = true)
end

function ODsolve_dense(F,x)
  @time EigenvalueSolver.solve_OD_dense(F,x)
end

function HCsolve(F)
  @time HomotopyContinuation.solve(F;threading=true)
end


function run()
  @polyvar x[1:3]

  f = [x[1]^2-x[2]^2,3x[1]^2-x[2]^5,2x[1]+9x[2]]

  @time HomotopyContinuation.solve(f;threading = true)
end

function convertgfptofloat(gfp)
  parse(Float64,string(gfg))
end


function solvegroebnerp(I,p,var)
  G = gens(I)
  sol = zeros(1:length(var))
  for i = 1:length(var) # maybe just up to length(var)
    duminput = zeros(Int,length(var))
    Gvars = vars(G[i])
    indices = Int[]

    for j = 1:length(Gvars)
    push!(indices,findfirst(isequal(Gvars[j]),var))
    end

    for k = 1:p
      duminput[indices[i]] = k
      if Singular.evaluate(G[i],duminput)== 0
        sol[indices[i]] = mod(k,5)
      else
        nothing
    end
  end
end
sol
end#solvegroeber

# evaluate function for an ideal
function evaluate(I::Oscar.MPolyIdeal, v)
    R = base_ring(I)
    return ideal([R(Oscar.evaluate(f,v)) for f in gens(I)])
end

# convert to univariate
function to_univar(a, var::Int, Kx)
    r = zero(Kx)
    for (c, ev) in zip(coefficients(a), exponent_vectors(a))
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
        elim = gens(eliminate(I, rest))[1]
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
        enum!(evaluate(I, vec), vcat(part, [v]), els, ans)
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
        enum!(evaluate(I, vec), part, els, ans)
    end
    ans
end

# example taken from the documentation of Magma
function example(char=7823)
    R,(x,y,z,w) = PolynomialRing(GF(char),["x","y","z","w"])
    I = ideal([4*x*z+2*x*w+y^2+4*y*w+7821*z^2+7820*w^2,4*x^2+4*x*y+7821*x*w+7822*y^2+7821*y*w+7821*z^2+7819*z*w+7820*w^2])
    @time length(proj_points(I))
end
# over an extension field
function exampleGF()
    R,(x,y) = PolynomialRing(FiniteField(19,3,"a")[1],["x","y"])
    I = ideal([y^27-x^8*(1-x)])
    @time length(points(I))
end




function netwons_method_qp(f,p)
  sol_modp = solve_mod(f,p)
  sol_mod_pn = sol_modp
  for n in 2:200
    sol_mod_pn = extend_mod(sol_modp,p,n)
  end
  final_sol = find_exact_sol(sol_mod_pn)
  double_check(f,final_sol)

  final_sol
end


function solve_mod(f,p,varinput)
  # use groebner basis from Oscar project
  # Nemo.jl, Singular.jl, or GroebnerBasis.jl

  # I have assumed input f is made through use of HomotopyContinuation

  # Do not delete the stuff below. It is incomplete
  # and does not do everything desired so far.
  # So far, one must define f using an array.
  # for instance, type in @polyvar x[1:3], f = [x[1]^2+1,x[2]+1,x[1]]
  # p=3, and varinput = x.
  # run the program. It should return 1 since the ideal generated by
  # f is the unit ideal.

  # for another test, try f=[x[1]^3-x[1]^2+2, x[2]-2x[1]^2+1, x[3]-3x[1]+5]
  # the output should be with generators
  # (x[3] - 1, x[1]^2 + x[2] + 1, x[1]*x[2] + x[1] - x[2], x[2]^2 - x[1] + 1)

  F = Fp(p)
  variablesnew = string.(varinput)
  R, (indeter) = Singular.PolynomialRing(F,variablesnew)
  #println(Singular.PolynomialRing(F,variablesnew))
  g = repeat([indeter[1]],length(f))
  for i in 1:length(f)
    g[i] = f[i]((varinput...,) => (indeter...,))
  end
  I = Ideal(R,g)
  I_groebner = std(I)
  # Need to solve mod p
end

function pseudo_inv_finite_field(M)
  # input M for matrix
  # formula from wikipedia: A^+ = (A^*A)^{-1}A^*
  # or A^*(A^*A)^{-1} depending on whether A has LI col/row res.

  if rank(M) == size(M)[1]
    pseudoinv = transpose(M)*inv(M*transpose(M))
  elseif rank(M) == size(M)[2]
    pseudoinv = inv(transpose(M) * M) * transpose(M)
  else
    pseudoinv = "cannot compute pseudoinverse of matrix input"
  end

  # only if our matrices are full rank is this well-defined
  # example: F[1 2; 2 3] gives 2x2 matrix with entries 1,2 ; 2,0.
  # Current issue is that M needs to be defined via above way before inputting.
  # this goes not compute M over modulo p. One must construct
  # M similarly to F[1 2; 2 4] by setting F = GF(p) for p prime.
  # Needs HomotopyContinuation to run.

  return pseudoinv
end

function jacobian_sol(f,sol,p)
  # makes Jacobian matrix of system f of polynomials and evalutes at sol
  # example: f =[x[1]^2+1,x[2]*x[1]]
  F = System(f)
  dummyjacob = jacobian(F,sol)
  for i in 1:ncols(dummyjacob)
    for j in 1:nrows(dummyjacob)
      dummyjacob[i,j] = dummyjacob[i,j] % p
    end
  end
  return dummyjacob
end

function evaluate(f,sol,xs)
  check(g) = double_check_single(g,sol,xs)
  # broadcasts the function to the array
  output = check.(f)
  return output
end

function extend_mod(sol,f,xs,m,p)
  # a_{n+1} = a_n - Df(a_n)^+ * f(a_n) Sketch of idea below
  # for stuff in so
  #   new_stuff(stuff)= stuff+ Df(stuff)^+ * f(stuff)
  #   sol = new_stuff.(stuff)
  # end
  # return sol
  # Need to have place in which this process terminates
  for n = 1:m
    extendsol = sol - pseudo_inv_finite_field(jacobian_sol(f,sol,p))*evaluate(f,sol,xs)
    sol = extendsol
  end
  return sol

  # Test to check this works. Run in repl the following.
  # @polyvar x[1:3], f = [x[1]+1, x[1]x[2]+2, x[1]+x[2]+x[3]],
  # sol = [1, 3, 5], test.extend_mod(sol,f,x,10,3).
  # observe sol does not satisfy f. however, the out put should be [-1, 2,-1]
  # which does satisfy f
end

function find_exact_sol(sol)
  # lindep (in Nemo)
end

function double_check_single(g,sol,xs)
  #println(g)
  # for n in 1:length(sol)
  #   g=subs(g, (vs...,)=>(sol...,))
  # end

  g( (xs...,) => (sol...,))
end

function double_check(f,sol,xs)
  # output = []
  # for g in f
  #   double_check_single(g,sol)
  #   push!(output, double_check_single(g,sol))
  # end
  # output

  # defines a function
  check(g) = double_check_single(g,sol,xs)
  # broadcasts the function to the array
  output = check.(f)
  # println(output)
  all(output .== zeros(length(output)))
end


end#module
