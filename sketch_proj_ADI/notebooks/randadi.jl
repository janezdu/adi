module RandAdi

# usings and imports
using Random
using LinearAlgebra
using SparseArrays
using Elliptic
# using Plots

# exported functions
export rand_matsolve, rand_linsolve, adi_solve, adi_parameters

function rand_linsolve(A,b,B, Stype, seed=nothing, verbose=false)
    """
    Solves using randomized kaczmarz.
    Randomly selects a row per iteration to update x with

    INPUT:  matrix A mxn,
            vector b mx1,
            matrix B nxn
            string Stype: "coordvec", ""

    OUTPUT: vector x nx1 that solves Ax = b
    """
    m,n = size(A)
    x0 = ones(n)
    sols = [x0] # list of x_0 ... x_n
    errs = zeros(0)
    if seed != nothing
        Random.seed!(seed)
    end
    @assert Stype == "coordvec"

    Binv = inv(B)
    xprev = x0

    while true
        i = rand(1:m)

        # S is a unit coordinate vector (temporary)
        S = zeros(m)
        S[i] = 1

        C = Binv * A' * S * pinv(S' * A * Binv * A' * S) * S'
        res = A * xprev - b
        xnew = xprev - C * res

        push!(sols, xnew)
        err = norm(A * xnew - b)
        push!(errs,err)
        if err < 1e-5
            break
        end

        xprev = xnew
    end

    if verbose
        return sols,errs
    else
        return sols[end]
    end
end

function rand_matsolve(A,B,seed=nothing)
    """
    Helper function that stacks together multiple Ax=b solvers
    INPUT:
            A: matrix mxn
            B: matrix mxk
    OUTPUT: X: matrix such that AX = B
    """
    
    m1,n = size(A)
    m2,k = size(B)
    sol = zeros(n,k)
    @assert m1==m2
    s(b) = rand_linsolve(A,b,I,"coordvec", seed)
    ret = mapslices(s,B,dims=[1])
    return ret
end

function adi_solve(A,B,F,N,p,q)
    """
    Uses randomized ADI to solve the equation AX - XB = F
    INPUT:
            A: matrix mxn
            B: matrix mxk
    OUTPUT: X: matrix such that AX = B
    """
    m,m2 = size(A)
    n,n2 = size(B)

    @assert m==m2
    @assert n==n2

    sols = []
    Xprev = zeros((m,n))
    for i = 1:N
        Ahalf = (A - p[i] * I)
        Bhalf = Xprev * (B - q[i] * I) + F
        Xhalf = rand_matsolve(Ahalf,Bhalf)

        Asolve = (B - q[i] * I)
        Bsolve = (A - q[i] * I) * Xhalf - F
        X = (rand_matsolve(Asolve', Bsolve'))'

        Xprev = X
        push!(sols,X)
    end
    return sols
end

function adi_parameters(a,b,c,d,J)
    p = zeros(J)
    q = zeros(J)
    gamma = abs((c-a) * (d-b)  / (c-b)*(d-a))
    alpha = -1 + 2*gamma + 2 * sqrt(gamma^2 - gamma)
    kappa = sqrt(1 - 1 / (alpha^2))

    TT(t) = (-alpha * t -1)/(t + alpha)

    for j = 1:J
        z = (2 * j + 1) / (2 * J) * Elliptic.K(kappa)
        t = -alpha * Elliptic.Jacobi.dn(z,kappa)
        p[j] = TT(t)
        q[j] = TT(-t)
    end

    return p,q
end


end
