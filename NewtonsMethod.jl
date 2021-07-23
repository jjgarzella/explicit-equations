module NewtonsMethod

using HomotopyContinuation
using Oscar
using Nemo
using Singular
using GroebnerBasis

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

  F = Fp(p)
  variablesnew = string.(varinput)
  R, (indeter) = Singular.PolynomialRing(F,variablesnew)
  #println(Singular.PolynomialRing(F,variablesnew))
  g = repeat([indeter[1]],length(f))
  for i in 1:length(f)
    g[i] = f[i]((varinput...,) => (indeter...,))
  end
  I = Ideal(R,g)
  I_groebner = f4(I)
  # Cannot get Singular.LibSolve.solve(I_groebner) to work... because it NEEDS characteristic 0
  # Singular.LibSolve.solve(I_groebner)
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
    pseudoinv = "cannot compute pseudoinverse of matrix input. Check rank."
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
  for n = 1:m
    extendsol = sol - pseudo_inv_finite_field(jacobian_sol(f,sol,p))*evaluate(f,sol,xs)
    sol = extendsol
  end
  return sol

  # Test to check this works. Run in repl the following.
  # @polyvar x[1:3], f = [x[1]+1, x[1]x[2]+2, x[1]+x[2]+x[3]],
  # sol = [1, 3, 5], NewtonsMethod.extend_mod(sol,f,x,10,3).
  # observe sol does not satisfy f. however, the out put should be [-1, 2,-1]
  # which does satisfy f
  # to check it works, run
  # NewtonsMethod.double_check(f,NewtonsMethod.extend_mod(sol,f,x,10,3),x)
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