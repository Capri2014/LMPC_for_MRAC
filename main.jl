# Below the packages need to run this code
using JuMP
using Ipopt
using JLD
using PyPlot

# Including some files
include("classes.jl")
include("LMPC_models.jl")
include("SolveLMPCProblem.jl")

# Delarins some types which are in the file classes.jl
SystemParams = TypeSystemParams()
LMPCparams   = TypeLMPCparams()
LMPCSol      = TypeLMPCSol()

# Now I define my dynamics
SystemParams.A  = [1.0 1.0; 0.0 1.0]
SystemParams.B  = [1.0 0.0; 1.0 0.0]
SystemParams.n  = 2
SystemParams.d  = 2

# Now define mpc parameters
LMPCparams.N = 4
LMPCparams.Q = [1.0 0.0; 0.0 1.0]
LMPCparams.R = [1.0 0.0; 0.0 1.0]

# This buils the optmization problem, basically Julia does
# some preprocessing to compile your problem
mdl    = LMPC_Model(LMPCparams,SystemParams)


# Initializing some variables
x_LMPC      = zeros(2,500)
x_LMPC[:,1] = [-1.0, -1.0]

u_LMPC      = zeros(2,500)
cost_LMPC   = ones(1,500)
t = 1

# This is the time loop
while cost_LMPC[t] > 0.000000003
    # At time t solve the optimization problem
    solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t])

    # Update your closed-loop 
    x_LMPC[:,t+1]  = LMPCSol.x[2,:]' # Note here no model missmatch
    u_LMPC[:,t]    = LMPCSol.u[1,:]'
    cost_LMPC[t+1] = LMPCSol.cost
    println("Cost of the finite time LQR problem at step ",t, " is ", cost_LMPC[t+1])

    t=t+1    
end

figure()
plot(x_LMPC[1,:]',x_LMPC[2,:]', "-ro")
    
grid(1)
title("Closed-loop Trajectory")
axis("equal")
xlabel("x_1")
ylabel("x_2")