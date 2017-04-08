using JuMP
using Ipopt
using JLD

include("classes.jl")
include("LMPC_models.jl")
include("SolveLMPCProblem.jl")

SystemParams = TypeSystemParams()
LMPCparams   = TypeLMPCparams()
LMPCSol      = TypeLMPCSol()


SystemParams.A  = [1.0 1.0; 0.0 1.0]
SystemParams.B  = [1.0 0.0; 1.0 0.0]
SystemParams.n  = 2
SystemParams.d  = 2

LMPCparams.N = 4
LMPCparams.Q = [1.0 0.0; 0.0 1.0]
LMPCparams.R = [1.0 0.0; 0.0 1.0]


mdl    = LMPC_Model(LMPCparams,SystemParams)



x_LMPC      = zeros(2,500)
x_LMPC[:,1] = [-1.0, -1.0]

u_LMPC      = zeros(2,500)
cost_LMPC   = ones(1,500)
t = 1

while cost_LMPC[t] > 0.000000003
    solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t])
    x_LMPC[:,t+1]  = LMPCSol.x[2,:]'
    u_LMPC[:,t]    = LMPCSol.u[1,:]'
    cost_LMPC[t+1] = LMPCSol.cost
    println("LMPC cost at step ",t, " is ", cost_LMPC[t+1])

    t=t+1    
end
