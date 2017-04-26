using JuMP
using Ipopt
using JLD
using PyPlot
using JLD
using PyPlot

include("classes.jl")
include("LMPC_models.jl")
include("SolveLMPCProblem.jl")
include("ComputeFeasibleTraj.jl")
include("ComputeCost.jl")


SystemParams = TypeSystemParams()
LMPCparams   = TypeLMPCparams()
LMPCSol      = TypeLMPCSol()

# Initialize System Parameters
n = 2
d = 2

Ar = [1.6 0.9; 
      0.0 1.8]

SystemParams.Ad  = [0.4 1.5; 
                    0.0 0.4]

SystemParams.Ar  = Ar

SystemParams.B  = [1.0 0.0; 
                   0.0 1.0]
SystemParams.n  = n
SystemParams.d  = d

LMPCparams.N  = 2
LMPCparams.Q  = [1.0 0.0; 
                 0.0 1.0]

LMPCparams.Qe = 1*[1.0 0.0; 
                   0.0 1.0]

LMPCparams.R  = [1.0 0.0; 
                 0.0 1.0]


x0          = [-38.0, -25.0]


# Compute First Feasible Iteration    
K_r    = -[0.9 0.7;
            0  1.6]


x_feasible, u_feasible = Feasible_Traj(SystemParams, x0, K_r)

# Initialize SS and Q function for first feasible iteration
Buffer = 200

SS   = zeros(n+n+3, Buffer+3, 20)
Qfun = zeros(1, Buffer+3, 20)
time = zeros(20)
time = round(Int64, time)

time[1] = size(x_feasible)[2] + 3 #This is to take into account the basis

x_LMPC              = zeros(n + n + 3,Buffer+3)
u_LMPC              = zeros(n + 3,Buffer-1+3)
x_real              = zeros(n , Buffer+3, 50)
u_real              = zeros(d,Buffer-1+3)
K_real        = zeros(2,2,Buffer+3)
K_real[:,:,1] = K_r


StateBasis = [0 0 0;
              0 0 0;
              0 0 0;
              0 0 0;
              1 0 0;
              0 1 0;
              0 0 1];

InputBasis = [0 0 0;
              0 0 0;
              0 0 0;
              0 0 0;
              0 0 0];

x_LMPC[:,1:time[1]]   = [x_feasible StateBasis] 
u_LMPC[:,1:time[1]-1] = [u_feasible InputBasis] 


it = 1
SS[:, 1:time[it], it]   = x_LMPC[:,1:time[it]]
Qfun[:, 1:time[it], it] = ComputeCost(x_LMPC[:,1:time[it]], u_LMPC[:,1:time[it]], LMPCparams)

# Now start with the Second iteration (The first is for the feasible trajectory)
it = 2
Difference = 1
while (abs(Difference) > (1e-7))&&(it<10)
    # Vectorize the SS and the Q function

    SSdim = sum(time)
    ConvSS   = zeros(2*n + 3, SSdim)
    ConvQfun = zeros(SSdim)
    Counter  = 1

    for ii = 1:it-1
        for kk = 1:time[ii]
            ConvSS[:,Counter]  = SS[:, kk, ii]
            ConvQfun[Counter]  = Qfun[1, kk, ii]

            Counter = Counter + 1
        end
    end

    # Here start the iteration
    x_LMPC      = zeros(n + n + 3,Buffer+3)
    x_LMPC[:,1] = x_feasible[:,1]

    K_real        = zeros(2,2,Buffer+3)
    K_real[:,:,1] = K_r

    x_real[:, 1, it]   = x_feasible[1:2,1]


    u_LMPC      = zeros(5,Buffer+3-1)
    cost_LMPC   = ones(1,500)

    # Define the model at the j-th iteration (Need to define it at each iterations as SS and Q function change)
    mdl    = LMPC_Model(LMPCparams,SystemParams, SSdim)
    
    # Enter the time loop for the LMPC at the j-th iteration
    t = 1
    while ((cost_LMPC[t] > (1e-5))&&(t<Buffer-1))
        
        if t == 1
            solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun) 
            solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun) 
        else
            solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun) 
            solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun) 
        end

        x_LMPC[:,t+1]  = LMPCSol.x[:,2]
        u_LMPC[:,t]    = LMPCSol.u[:,1]

        K_real[:,:,t+1]  = K_real[:,:,t] + [u_LMPC[3,t] u_LMPC[4,t]; 0 u_LMPC[5,t]]

        u_real[:,t]       = K_real[:,:,t+1]*x_real[:,t, it] + u_LMPC[1:2,t]
        x_real[:,t+1,it]  = Ar * x_real[:,t, it] + u_real[:,t]  

        x_LMPC[3:4,t+1] = x_real[:,t+1,it] - x_LMPC[1:2,t+1]
        
        cost_LMPC[t+1] = LMPCSol.cost
        println("LMPC cost at step ",t, " of iteration ", it," is ", cost_LMPC[t+1])

        t=t+1    
    end

    # Now post process the data after the LMPC has converged
    time[it] = t + 3 # This is to take into account the basis
    x_LMPC[:,t+1:t+3] = StateBasis
    u_LMPC[:,t+1:t+3] = InputBasis

    # Add data to SS and Q function
    SS[:, 1:time[it], it]   = x_LMPC[:,1:time[it]]
    Qfun[:, 1:time[it], it] = ComputeCost(x_LMPC[:,1:time[it]], u_LMPC[:,1:time[it]], LMPCparams)

    Difference = Qfun[1,1,it-1]-Qfun[1,1,it]

    it = it + 1

end
it = it - 1
for i = 1:it
    println(i,"-th itearion cost; ", Qfun[1,1,i])
end


for i = 2:it
    figure()
    plot(x_feasible[1,:]',x_feasible[2,:]', "-ro")
    hold(1)

    j = 1
    plot(SS[1, 1:time[j], j]', SS[2, 1:time[j], j]', "-go", label="Safe Set")
    for j = 2:(i-1)
        plot(SS[1, 1:time[j], j]', SS[2, 1:time[j], j]', "-go")
    end

    plot(SS[1, 1:time[i], i]', SS[2, 1:time[i], i]', "-ko", label="Trajectory Reference System")
    plot(x_real[1, 1:time[i], i]', x_real[2, 1:time[i], i]', "-r*", label="Trajectory Real System" )
    
    grid(1)
    LMPCIteration = i -1
    title("Closed-loop trajectory at iteration $LMPCIteration")
    axis("equal")
    xlabel(L"$x_1$", size=24)
    ylabel(L"$x_2$", size=24)
    legend(loc="lower right",fancybox="true")
end

figure()
hold(1)
i = 2
vec1 = SS[3, 1:time[i], i].^2 + SS[4, 1:time[i], i].^2
println("Here ",vec)
plot(1:time[i], vec1[:] , "-ro", label="1st Iteration")

i = it
vecss = SS[3, 1:time[i], i].^2 + SS[4, 1:time[i], i].^2
println("Here ",vec)
plot(1:time[i], vecss[:] , "-g*", label="Steady State")

grid(1)
title("LMPC Steady State")
axis("equal")
xlabel(L"Time step", size=24)
ylabel(L"$||x_2||_2^2$", size=24)
legend()
