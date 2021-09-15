module editrels
# this was written to parse through the file separatedBFrels/Rels23withrr
# from the paper New Explicit Constructions of Surfaces of General Type
# and prepare for it to be inputted into the HomotopyContinuation.jl solver

using HomotopyContinuation
using DynamicPolynomials

# this sets up the file by taking in the file with equations
# and then outputting a new file which can then be copy pasted into the
# repl and applied to the solve function of HomotopyContinuation
function setup(filename,outputfilename)
    open(filename,"r") do io
    open(outputfilename,"w") do jo
        for line in eachline(io)
            line = replace(line, "{" => "[")
            line = replace(line, "}" => "]")
            for i =1:6
                line = replace(line, "d$i" => "d[$i]")
            end
            for j=10:19
                line = replace(line, "c$j" =>"c[$(j+1)]")
            end
            # Need to do these replacements first otherwise it messes things up
            for j=0:9
                line = replace(line, "c$j" =>"c[$(j+1)]")
            end
            line = replace(line, "rr" => "sqrt(complex(-15))")
            write(jo,line)
        end
    end;
    end;
end#function setup(filename)

"""
setup("separatedBFrels/Rels23withrrIn.txt","separatedBFrels/Rels23withrrOut.txt")
"""

function setupmod(filename,outputfilename)
    open(filename,"r") do io
    open(outputfilename,"w") do jo
        for line in eachline(io)
            line = replace(line, "{" => "[")
            line = replace(line, "}" => "]")
            for i =1:6
                line = replace(line, "d$i" => "d[$i]")
            end
            for j=10:19
                line = replace(line, "c$j" =>"c[$(j+1)]")
            end
            # Need to do these replacements first otherwise it messes things up
            for j=0:9
                line = replace(line, "c$j" =>"c[$(j+1)]")
            end
            line = replace(line, "rr" => "2")
            write(jo,line)
        end
    end;
    end;
end#function setup(filename)

"""
editrels.setupmod("separatedBFrels/Rels23withrrIn.txt","separatedBFrels/Rels23withrrOut.txt")
"""


function simplify(filename,outputfilename)
    open(filename,"r") do io
    open(outputfilename,"w") do jo
    for line in eachline(io)
        line = replace(line, ", " => "]
["
        #line = replace(line, ", " => "]
#A = [")
        write(jo,line)
    end
    end;
    end;
end

"""
editrels.simplify("separatedBFrels/Rels23withrrOut.txt","separatedBFrels/Rels23Simp.txt")
"""

function reindex(filename,outputfilename)
    open(filename,"r") do io
    open(outputfilename,"w") do jo
    n = 0
    for line in eachline(io)
        n += 1
        line = replace(line, "A" => "A$n")
        write(jo,line)
    end
    end;
    end;
end

"""
editrels.reindex("separatedBFrels/Rels23Simp.txt","separatedBFrels/Rels23Final.txt")
"""


function reindex2(filename,outputfilename)
    @polyvar c[1:20]
    @polyvar d[1:6]
    open(filename,"r") do io
    open(outputfilename,"w") do jo
    n = 0
    for line in eachline(io)
        n += 1
        #line = replace(line, "A" => "
#A$n")
        write(jo,line)
    end
    end;
    end;
end


function separatefiles(filename)
    n = 0
    open(filename,"r") do io
    for line in eachline(io)
        n += 1
    open("separatedBFrels/RelsEq$n","w") do jo
        write(jo,line)
    end
    end;
    end;
end

function separatefilesmod(filename)
    n = 0
    open(filename,"r") do io
    for line in eachline(io)
        n += 1
    open("separatedBFrels/RelsEqMod$n","w") do jo
        write(jo,line)
    end
    end;
    end;
end


"""
editrels.reindex2("separatedBFrels/Rels23Simp.txt","separatedBFrels/Rels23Final.txt")
Edit the first entry
editrels.separatefiles("separatedBFrels/Rels23Final.txt")

function read(s::IO, ::Type{Complex{T}}) where T<:Real
    r = read(s,T)
    i = read(s,T)
    Complex{T}(r,i)
end
"""

# this function will merge the equations
# in each file together

function convdynam(poly)
    #convert(Vector{DynamicPolynomials.Polynomial{true, ComplexF64}}, poly)
end

function mergerels()
    @polyvar c[1:20]
    @polyvar d[1:6]
    S = []
    for n = 1
        #S = open(f->read(f, Vector{Polynomial{true, ComplexF64}} ),"separatedBFrels/RelsEq$n")
        open(f->vcat(S,eval(Meta.parse(f))),"separatedBFrels/RelsEq$n")
    end
    println(S)
end

"""
editrels.mergerels()
convert(Vector{DynamicPolynomials{true, ComplexF64}}, input)
"""

# function setup(filename)
#     open(filename,"r") do io
#     open("output.txt","w") do jo
#         for line in eachline(io)
#             for i =1:6
#                 replace(line, "d$i" => "d[$i]")
#             end
#             for j=10:19
#                 replace(line, "c$j" =>"c[$(j+1)]")
#             end
#             for j=0:9
#                 replace(line, "c$j" =>"c[$(j+1)]")
#             end
#         write(jo,line)
#         end
#     end;
#     end;
# end#function setup(filename)

# function setup(rangeofc,rangeofd)
#     @polyvar c[1:rangeofc]
#     @polyvar d[1:rangeofd]
#     cdum = fill("s",rangeofc)
#     ddum = fill("s",rangeofd)
#     for i = 1:rangeofc
#         cdum[i] = ("c"*string(i))
#     end
#     for i = 1:rangeofc
#         c[i]=cdum[i]
#     end
#     for i = 1:rangeofd
#         ddum[i] = string.("d"*string(i))
#     end
#     for i = 1:rangeofd
#         d[i]=ddum[i]
#     end
#     rr = sqrt(complex(-15))
# end
#
# function setup
#     (tmppath, tmpio) = mktemp()
#     open(file) do io
#         for line in eachline(io, keep=true) # keep so the new line isn't chomped
#             if predicate(line)
#                 line = modifier(line)
#             end
#             write(tmpio, line)
#         end
#     end
#     close(tmpio)
#     mv(tmppath, file, force=true)
# end



end #module
