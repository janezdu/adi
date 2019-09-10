include("../src/adi.jl")
using .Adi, Test

#--- Small Test
A = rand(Float64,5,5)
B = rand(Float64,3,3)
F = rand(Float64,5,3)

N = 5
p = ones(N)
q = ones(N)

sols = Adi.adi_solve(A,B,F,N,p,q)


#---


#--- Medium Random Test
A = rand(Float64,25,25)
B = rand(Float64,15,15)
F = rand(Float64,25,15)

N = 5
p = ones(N)
q = ones(N)

sols = adi(A,B,F,N,p,q)

display("text/plain",sols[end])
# X = sols[end]
# println(norm(A *X - X * B))
