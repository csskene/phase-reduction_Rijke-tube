using SparseArrays
using LinearAlgebra
using QuadGK

mutable struct modelParameters
    N
    c1
    c2
    K
    xf
    τ
    A1
    solDir
    T
    tMax
    solAdj
    phaseSensNorm
end

function setParams(N,c1,c2,K,xf,τ,print=false)
	D  = zeros(N,N)
	Om = zeros(N,N)
	Id = zeros(N,N)
    for i in 1:N
        D[i,i]  = 2*damping(c1,c2,i)*i*π
        Om[i,i] = i*i*π*π
        Id[i,i] = 1
    end
    A1 = [zeros(N,N) Id;-Om -D]
    A1 = sparse(A1)
    params = modelParameters(N,c1,c2,K,xf,τ,A1,false,false,false,false,false)
    if(print)
	println("#####################################")
	println("## Rijke tube Galerkin modes model ##")
	println("#####################################")
	println("\n##  Parameters ##")
	println("N  = ",params.N)
	println("c1 = ",params.c1)
	println("c2 = ",params.c2)
	println("K  = ",params.K)
	println("xf = ",params.xf)
	println("τ  = ",params.τ)
	println("#################\n")
    end
    return params
end

function damping(c1,c2,i)
    return (1/(2π))*(c1*i+c2*sqrt(1/i))
end

function flameVel(xd,p)
    uf = 0
    for i in 1:p.N
        uf = uf + xd[i]*cos(i*π*p.xf)
    end
    return uf
end

function RHS(u,h,p,t)

    uf = flameVel(h(p,t-p.τ),p)

    y = p.A1*u

    for i in 1:p.N
        y[N+i] = y[N+i] - i*π*p.K*(sqrt(abs(1/3+uf))-sqrt(1/3))*sin(i*π*p.xf)
    end
    return y
end

function pressureVelocity(modes,x,params)
    # Returns the pressure and velocity from the Galerkin modes at x
    u = zeros(length(x),)
    p = zeros(length(x),)
    for i in 1:params.N
	    u += modes[i]*cos.(i*π*x)
        p -= modes[i+params.N]*sin.(i*π*x)*γ*Ma/(i*π)
    end

    return u,p
end

function energy(modes,params)
    # Calculates the non-dimensional acoustic energy per unit volume
    # from the Galerkin modes
    E = 0
    for i=1:params.N
        E += modes[i]^2+(modes[params.N+i]/(i*π))^2
    end
    E *= 0.5
    return E
end


## Linear and Adjoint routines
function linHeatRelease(t,p)
    # Creates the time-delayed system matrix A2 from the direct solution at t.
    Q = zeros(p.N,p.N)
    z = p.solDir(t-p.τ)
    ufb = flameVel(z,p)

    for i in 1:p.N
        for j in 1:p.N
            Q[i,j]=0.5*i*π*p.K*sin(i*π*p.xf)*sign(1/3+ufb)*cos(j*π*p.xf)/sqrt(abs(1/3+ufb))
        end
    end
    A2 = -[zeros(p.N,p.N) zeros(p.N,p.N);Q zeros(p.N,p.N)]
    return A2
end

function linHeatReleaseMat(base,p)
    # Creates the linearised matrix A2 from the state base
    Q = zeros(p.N,p.N)
    z = base
    ufb = flameVel(z,p)

    for i in 1:p.N
        for j in 1:p.N
            Q[i,j]=0.5*i*π*p.K*sin(i*π*p.xf)*sign(1/3+ufb)*cos(j*π*p.xf)/sqrt(abs(1/3+ufb))
        end
    end
    A2 = -[zeros(p.N,p.N) zeros(p.N,p.N);Q zeros(p.N,p.N)]
    return A2
end

function adjTime(tAdj)
    # Allows us to solve the adjoint forwards in time
    return params.tMax .- tAdj
end

function RHSadj(u,h,p,t)
    A2 = linHeatRelease(adjTime(t-p.τ),p)
    return p.A1'*u + A2'*h(p,t-p.τ)
end

# Misc functions

function bilinearForm(t,p)
    # Inner product part
    innerProduct = dot(p.solAdj(adjTime(t)),p.solDir(t,Val{1}))
    # Integral part
    integrand(ξ) = dot(p.solAdj(adjTime(t+ξ+p.τ)),linHeatRelease(t+ξ+p.τ,p)*p.solDir(t+ξ,Val{1}))
    integral,error = quadgk(integrand, -p.τ, 0., rtol=1e-4)

    norm = innerProduct + integral
    return norm, innerProduct, integral
end

function phaseCouplingFunction(ϕ,ωf,f,p)
    Tf = 2*π/ωf
    integrand(s) = dot(p.solAdj(adjTime(ϕ+params.T*s/Tf)),f(s))
    integral,error = quadgk(integrand, 0, Tf, rtol=1e-4)
    Γ = integral / Tf
    return Γ
end

function phaseCouplingFunctionmn(ϕ,ωf,f,m,n,p)
    Tf = 2*π/ωf
    integrand(s) = dot(p.solAdj(adjTime(ϕ+(n/m)*params.T*s/Tf)),f(s))
    integral,error = quadgk(integrand, 0, m*Tf, rtol=1e-4)
    Γ = integral / (m*Tf)
    return Γ
end
