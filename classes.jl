type TypeSystemParams          # parameters for MPC solver
    g::Float64
    rho::Float64
    dt::Float64
    xF::Array{Float64,2}
    TypeSystemParams( g = 9, rho = 4,  dt= 0, xF = [0 0 0 0] ) = new(g, rho, dt, xF)
end

type TypeLMPCparams          # parameters for MPC solver
    Q::Array{Float64,2}
    Qe::Array{Float64,2}
    R::Array{Float64,2}
    N::Int64
    TypeLMPCparams( Q=[1 1;1 1], Qe=[1 1;1 1], R=[1 1;1 1], N = 0 ) = new(Q, Qe, R, N)
end

type TypeLMPCSol          # parameters for MPC solver
    x::Array{Float64,2}
    u::Array{Float64,2}
    cost::Float64
    TypeLMPCSol( x=[1 1;1 1], u=[1 1;1 1], cost = 0) = new(x, u, cost)
end
