include("../src/adi.jl")
using .Adi, Test, LinearAlgebra

#--- Small Test
Q,R = qr(rand(5,5))
Aevals = [1., 2., 3., 4., 5.]
Bevals = [11., 12., 13., 14., 15.]
A = Q * Diagonal(Aevals) * Q.T
B = Q * Diagonal(Bevals) * Q.T

F = rand(Float64,5,5)

a = minimum(Aevals)
b = maximum(Aevals)
c = minimum(Bevals)
d = maximum(Bevals)
N = 5
p,q = adi_parameters(a,b,c,d,N)

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
