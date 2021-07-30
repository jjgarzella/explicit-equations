module NewtonsMethod

using HomotopyContinuation
using Oscar
using Nemo
using Singular
using GroebnerBasis

function netwons_method_qp(f,p,var)
  sol_modp = solve_mod(f,p,var)
  sol_mod_pn = sol_modp
  for n in 2:200
    sol_mod_pn = extend_mod(sol_modp,p,n)
  end
  final_sol = find_exact_sol(sol_mod_pn)
  double_check(f,final_sol)

  final_sol
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

function solve_mod(f,p,var)
  # Example, type in @polyvar x[1:3], f = [x[1]^2+1,x[2]+1,x[1]]
  # p=3, and var = x.
  # the program should fail since this system cannot be solved

  # for another test, try f=[x[1]^3-x[1]^2+2, x[2]-2x[1]^2+1, x[3]-3x[1]+5]
  # the Groebner basis is (x[3] - 1, x[1]^2 + x[2] + 1, x[1]*x[2] + x[1] - x[2], x[2]^2 - x[1] + 1)
  # and the program fails to solve this

  # for a currently successful test, try f = [x[1]^3-2x[1]*x[2],x[1]^2*x[2]-2x[2]^2+x[1]]
  # which gives Groebner basis (x[1]^2, x[2]^2 + 2*x[1], x[1]*x[2])
  # and solution [0.0,0.0] in type Float64.

  F = Fp(p)
  variablesnew = string.(var)
  R, (indeter) = Singular.PolynomialRing(F,variablesnew)
  #println(Singular.PolynomialRing(F,variablesnew))
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

function pseudo_inv_finite_field(M)
  # input M for matrix
  # formula from wikipedia: A^+ = (A^*A)^{-1}A^*
  # or A^*(A^*A)^{-1} depending on whether A has LI col/row res.

  if rank(M) == size(M)[1]
    pseudoinv = transpose(M)*inv(M*transpose(M))
  elseif rank(M) == size(M)[2]
    pseudoinv = inv(transpose(M) * M) * transpose(M)
  else
    println("Cannot compute pseudoinverse of matrix input. Check rank.")
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
  for i in 1:nrows(dummyjacob)
    for j in 1:ncols(dummyjacob)
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
    extendsol = sol - pseudo_inv_finite_field(jacobian_sol(f,sol,p))*NewtonsMethod.evaluate(f,sol,xs)
    sol = extendsol
  end
  # extendsol does not do its work over mod p
  # observe in the example below that the output [-1,2,-1] is actually a Float64
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
