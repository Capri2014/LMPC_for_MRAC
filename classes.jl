type TypeSystemParams          # parameters for MPC solver
    Ar::Array{Float64,2}
    Ad::Array{Float64,2}
    B::Array{Float64,2}
    n::Int64
    d::Int64
    TypeSystemParams( Ar=[1 1;1 1], Ad=[1 1;1 1], B=[1 1;1 1], n = 0, d= 0 ) = new(Ar, Ad, B, n, d)
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
