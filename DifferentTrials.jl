function DifferentTrials()


include("classes.jl")
include("LMPC_models.jl")
include("SolveLMPCProblem.jl")
include("NominalLMPC_models.jl")
include("SolveNominalLMPCProblem.jl")
include("ComputeFeasibleTraj.jl")
include("ComputeCost.jl")
include("SystemID_Inloop.jl")
include("SystemID_Outloop.jl")

SystemParams = TypeSystemParams()
LMPCparams   = TypeLMPCparams()
LMPCSol      = TypeLMPCSol()
NLMPCSol     = TypeLMPCSol()

# Initialize System Parameters
close("all")
n = 2
d = 1

Ar = [1.6 0.9; 
      0.0 1.8]

Ar = [0.6 0.9; 
      0.0 0.8]

# Compute First Feasible Iteration    
# K_r    = -[1.1396; 2.2270]
# K_r    = -[2.0679; 2.9627]
K_r    = -[0.5; 0.7]

B  = [0.0, 1.0]

SystemParams.Ar  = Ar

SystemParams.B  = B
SystemParams.n  = n
SystemParams.d  = d


LMPCparams.N  = 4
LMPCparams.Q  = [1.0 0.0; 
                 0.0 1.0]

LMPCparams.Qe = 1*[1.0 0.0; 
                   0.0 1.0]

LMPCparams.R  = [1.0 0.0; 
                 0.0 1.0]


x0          = [-38.0, -25.0]




x_feasible, u_feasible, MeanEstimate, MSE = Feasible_Traj(SystemParams, x0, K_r)
x_feasible, u_feasible, NMeanEstimate, NMSE = Feasible_Traj(SystemParams, x0, K_r)

# Initialize SS and Q function for first feasible iteration
Buffer = 200

OptX = zeros(n, Buffer, 20)
OptU = zeros(1, Buffer, 20)

SS   = zeros(n, Buffer, 20)
NSS   = zeros(n, Buffer, 20)
OldU = zeros(1, Buffer, 20)
OldNU = zeros(1, Buffer, 20)
Qfun = zeros(1, Buffer, 20)
NQfun = zeros(1, Buffer, 20)
OptQ = zeros(1, Buffer, 20)
time = zeros(20)
time = round(Int64, time)

SaveMean = zeros(2,3,20)
SaveVari = zeros(1,2,20)
NSaveMean = zeros(2,3,20)
NSaveVari = zeros(1,2,20)

time[1] = size(x_feasible)[2]

figure()
plot(x_feasible[1,:]',x_feasible[2,:]', "-ro")

x_LMPC              = zeros(2,Buffer)
u_LMPC              = zeros(1,Buffer)
x_NLMPC              = zeros(2,Buffer)
u_NLMPC              = zeros(1,Buffer)

x_LMPC[:,1:time[1]]   = x_feasible
u_LMPC[:,1:time[1]-1] = u_feasible 
x_NLMPC[:,1:time[1]]   = x_feasible
u_NLMPC[:,1:time[1]-1] = u_feasible 

OptX[:,1:time[1]]     = x_feasible
OptU[:,1:time[1]-1]   = u_feasible 


it = 1
SS[:, 1:time[it], it]   = x_LMPC[:,1:time[it]]
NSS[:, 1:time[it], it]   = x_NLMPC[:,1:time[it]]

OldU[:, 1:time[it]-1, it]   = u_LMPC[:,1:time[it]-1]
OldNU[:, 1:time[it]-1, it]   = u_NLMPC[:,1:time[it]-1]
Qfun[:, 1:time[it], it]  = ComputeCost(x_LMPC[:,1:time[it]], u_LMPC[:,1:time[it]], LMPCparams)
NQfun[:, 1:time[it], it] = ComputeCost(x_NLMPC[:,1:time[it]], u_NLMPC[:,1:time[it]], LMPCparams)
OptQ[:, 1:time[it], it] = ComputeCost(x_LMPC[:,1:time[it]], u_LMPC[:,1:time[it]], LMPCparams)

# Now start with the Second iteration (The first is for the feasible trajectory)
it = 2
Difference = 1
while (abs(Difference) > (1e-7))&&(it<10)
    # Vectorize the SS and the Q function
    SSdim = sum(time)
    ConvSS   = zeros(n, SSdim)
    ConvQfun = zeros(SSdim)
    Counter  = 1

    for ii = 1:it-1
        for kk = 1:time[ii]
            ConvSS[:,Counter]  = SS[:, kk, ii]
            ConvQfun[Counter]  = Qfun[1, kk, ii]

            Counter = Counter + 1
        end
    end
    
    NSSdim = sum(time)
    NConvSS   = zeros(n, NSSdim)
    NConvQfun = zeros(NSSdim)
    Counter  = 1

    for ii = 1:it-1
        for kk = 1:time[ii]
            NConvSS[:,Counter]  = NSS[:, kk, ii]
            NConvQfun[Counter]  = NQfun[1, kk, ii]

            Counter = Counter + 1
        end
    end

    # Here start the iteration
    x_LMPC      = zeros(n,Buffer)
    x_LMPC[:,1] = x_feasible[:,1]
    x_NLMPC      = zeros(n,Buffer)
    x_NLMPC[:,1] = x_feasible[:,1]

    OptX[:,1,it]= x_feasible[:,1]

    u_LMPC      = zeros(1,Buffer)
    cost_LMPC   = ones(1,500)

    # Define the model at the j-th iteration (Need to define it at each iterations as SS and Q function change)
    mdl    =        LMPC_Model(LMPCparams, SystemParams,  SSdim,  MeanEstimate,  MSE)
    Nmdl   = NominalLMPC_Model(LMPCparams, SystemParams, NSSdim, NMeanEstimate)
    
    # ========================================================================================================
    # ===================== Enter the time loop for the LMPC at the j-th iteration ===========================
    # ========================================================================================================
    t = 1
    Max_x = 100000

    # println("MeanEstimate", MeanEstimate, "Variance ", MSE)
    # println("Nominal MeanEstimate", NMeanEstimate, "Nominal Variance ", NMSE)

    SaveMean[:,:,it] = MeanEstimate
    SaveVari[:,:,it] = MSE
    NSaveMean[:,:,it] = NMeanEstimate
    NSaveVari[:,:,it] = NMSE
    
    while t<20#((Max_x > (10))&&(t<Buffer-1))
        # 
        if t == 1
            #solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun, MeanEstimate, MSE) 
            solveLMPCProblem( mdl, LMPCSol,  x_LMPC[:,t],  ConvSS,  ConvQfun,  MeanEstimate,  MSE) 
            solveNominalLMPCProblem(Nmdl,NLMPCSol, x_NLMPC[:,t], NConvSS, NConvQfun, NMeanEstimate) 
        else
            #solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun, MeanEstimate, MSE) 
            # println(MeanEstimate, MSE)
            # println(NMeanEstimate, NMSE)
            
            solveLMPCProblem( mdl, LMPCSol,  x_LMPC[:,t],  ConvSS,  ConvQfun,  MeanEstimate,  MSE) 
            solveNominalLMPCProblem(Nmdl,NLMPCSol, x_NLMPC[:,t], NConvSS, NConvQfun, NMeanEstimate) 
        end
        Noise = 2*[2*randn(), 3*randn()]
        
        u_LMPC[:,t]   =  LMPCSol.u[:,1]
        u_NLMPC[:,t]  = NLMPCSol.u[:,1]

        x_LMPC[:,t+1]   = Ar * x_LMPC[:,t]  + *([0;1], u_LMPC[1,t]) + Noise#[Noise[1]*x_LMPC[1,t]; Noise[2]*x_LMPC[2,t]]*0.1
        x_NLMPC[:,t+1]  = Ar * x_NLMPC[:,t] + *([0;1], u_NLMPC[1,t]) + Noise#[Noise[1]*x_NLMPC[1,t]; Noise[2]*x_NLMPC[2,t]]*0.1

        OptU[1,t,it] = - dot([1.49455, 2.50961], OptX[:,t,it])
        OptX[:,t+1,it]  = Ar * OptX[:,t,it] + [0;1]* OptU[1,t,it] + Noise#[Noise[1]*OptX[1,t,it]; Noise[2]*OptX[2,t,it]]*0.1

        Max_x = max(abs(x_LMPC[1,t+1]), abs(x_LMPC[2,t+1]) )
        cost_LMPC[t+1] = LMPCSol.cost
        # println("LMPC cost at step ",t, " of iteration ", it," is ", cost_LMPC[t+1], " and Estimte is ", LMPCSol.a)
        # println("Value ", x_LMPC[:,t+1])

        # System ID at time t
        MeanEstimate,   MSE = SystemID_Inloop( time,t, it, SS,  x_LMPC,  OldU,  u_LMPC)
        NMeanEstimate, NMSE = SystemID_Inloop( time,t, it, NSS,  x_NLMPC,OldNU, u_NLMPC)

        # println("Estimte is ", MeanEstimate)
        t=t+1
    end
    # ========================================================================================================
    # ==================================== Back to the iterations loop =======================================
    # ========================================================================================================

    # Now post process the data after the LMPC has converged
    time[it] = t
    # Add data to SS and Q function
    OldU[:, 1:time[it]-1, it] = u_LMPC[:,1:time[it]-1]
    OldNU[:, 1:time[it]-1, it] = u_NLMPC[:,1:time[it]-1]
    SS[:, 1:time[it], it]     = x_LMPC[:,1:time[it]]
    NSS[:, 1:time[it], it]    = x_NLMPC[:,1:time[it]]

    NQfun[:, 1:time[it], it] = ComputeCost(x_NLMPC[:,1:time[it]], u_NLMPC[:,1:time[it]], LMPCparams)
    Qfun[:, 1:time[it], it] = ComputeCost(x_LMPC[:,1:time[it]], u_LMPC[:,1:time[it]], LMPCparams)
    OptQ[:, 1:time[it], it] = ComputeCost(OptX[:,1:time[it],it], OptU[:,1:time[it],it], LMPCparams)



    Difference = Qfun[1,1,it-1]-Qfun[1,1,it]

    MeanEstimate, MSE = SystemID_Outloop(time, it, SS, OldU)

    NMeanEstimate, NMSE = SystemID_Outloop(time, it, NSS, OldNU)
    it = it + 1

end
it = it - 1

# figure()
# hold(1)
# for i = 1:it-1
#     plot(SS[1,:,i]',SS[2,:,i]', "-ks", label="SS")
# end
# plot(x_LMPC[1,:]',x_LMPC[2,:]', "-ro", label="x_LMPC")
# plot(OptX[1,:,it]',OptX[2,:,it]', "-g*", label="Opt")
# legend()

figure()
index = it
plot(SS[1,:,index]',SS[2,:,index]', "-ko", label="x_LMPC at last it")
plot(OptX[1,:,index]',OptX[2,:,index]', "-g*", label="Opt")
legend()

figure()
index = 2
plot(SS[1,:,index]',SS[2,:,index]', "-ko", label="x_LMPC it 2")
plot(NSS[1,:,index]',NSS[2,:,index]', "-ro", label="x_NLMPC it 2")
plot(OptX[1,:,index]',OptX[2,:,index]', "-g*", label="Opt")
legend()


for i = 2:it-1
    println(i,"-th itearion cost LMPC; ", Qfun[1,1,i], " Nominal ", NQfun[1,1,i], "LMPC-Nominal ", Qfun[1,1,i]-NQfun[1,1,i])
    println(i,"-th itearion cost; ", OptQ[1,1,i])
end
Val = zeros(20,1)
for i = 2:it-1
    Val[i,1] = Qfun[1,1,i]-NQfun[1,1,i]
end

    return Val
end