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
d = 1

Ar = [1.6 0.9; 
      0.0 1.8]

B  = [0.0, 1.0]

SystemParams.Ar  = Ar

SystemParams.B  = B
SystemParams.n  = n
SystemParams.d  = d


LMPCparams.N  = 3
LMPCparams.Q  = [1.0 0.0; 
                 0.0 1.0]

LMPCparams.Qe = 1*[1.0 0.0; 
                   0.0 1.0]

LMPCparams.R  = [1.0 0.0; 
                 0.0 1.0]


x0          = [-38.0, -25.0]


# Compute First Feasible Iteration    
K_r    = -[1.1396; 2.2270]

x_feasible, u_feasible, MeanEstimate, MSE = Feasible_Traj(SystemParams, x0, K_r)

# Initialize SS and Q function for first feasible iteration
Buffer = 200

OptX = zeros(n, Buffer, 20)
OptU = zeros(1, Buffer, 20)

SS   = zeros(n, Buffer, 20)
OldU = zeros(1, Buffer, 20)
Qfun = zeros(1, Buffer, 20)
OptQ = zeros(1, Buffer, 20)
time = zeros(20)
time = round(Int64, time)

SaveMean = zeros(2,2,20)
SaveVari = zeros(1,2,20)

time[1] = size(x_feasible)[2]

figure()
plot(x_feasible[1,:]',x_feasible[2,:]', "-ro")

x_LMPC              = zeros(2,Buffer)
u_LMPC              = zeros(1,Buffer)

x_LMPC[:,1:time[1]]   = x_feasible
u_LMPC[:,1:time[1]-1] = u_feasible 
OptX[:,1:time[1]]     = x_feasible
OptU[:,1:time[1]-1]   = u_feasible 


it = 1
SS[:, 1:time[it], it]   = x_LMPC[:,1:time[it]]
OldU[:, 1:time[it]-1, it]   = u_LMPC[:,1:time[it]-1]
Qfun[:, 1:time[it], it] = ComputeCost(x_LMPC[:,1:time[it]], u_LMPC[:,1:time[it]], LMPCparams)
OptQ[:, 1:time[it], it] = ComputeCost(x_LMPC[:,1:time[it]], u_LMPC[:,1:time[it]], LMPCparams)

# Now start with the Second iteration (The first is for the feasible trajectory)
it = 2
Difference = 1
while (abs(Difference) > (1e-7))&&(it<20)
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

    # Here start the iteration
    x_LMPC      = zeros(n,Buffer)
    x_LMPC[:,1] = x_feasible[:,1]

    OptX[:,1,it]= x_feasible[:,1]

    u_LMPC      = zeros(1,Buffer)
    cost_LMPC   = ones(1,500)

    # Define the model at the j-th iteration (Need to define it at each iterations as SS and Q function change)
    mdl    = LMPC_Model(LMPCparams,SystemParams, SSdim, MeanEstimate, MSE)
    
    # ========================================================================================================
    # ===================== Enter the time loop for the LMPC at the j-th iteration ===========================
    # ========================================================================================================
    t = 1
    Max_x = 100000

    println("MeanEstimate", MeanEstimate, "Variance ", MSE)

    SaveMean[:,:,it] = MeanEstimate
    SaveVari[:,:,it] = MSE

    while t<40#((Max_x > (10))&&(t<Buffer-1))
        # 
        if t == 1
            #solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun, MeanEstimate, MSE) 
            solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun, MeanEstimate, MSE) 
        else
            #solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun, MeanEstimate, MSE) 
            solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun, MeanEstimate, MSE) 
        end
        Noise = 1*[2*randn(); 3*randn()]
        
        u_LMPC[:,t]    = LMPCSol.u[:,1]
        x_LMPC[:,t+1]  = Ar * x_LMPC[:,t] + *([0;1], u_LMPC[1,t]) + Noise

        OptU[1,t,it] = - dot([1.49455, 2.50961], OptX[:,t,it])
        OptX[:,t+1,it]  = Ar * OptX[:,t,it] + [0;1]* OptU[1,t,it] + Noise

        Max_x = max(abs(x_LMPC[1,t+1]), abs(x_LMPC[2,t+1]) )
        cost_LMPC[t+1] = LMPCSol.cost
        println("LMPC cost at step ",t, " of iteration ", it," is ", cost_LMPC[t+1], " and Estimte is ", LMPCSol.a)
        println("Value ", x_LMPC[:,t+1])

        # System ID
        SSdim_ID = sum(time) + t
        vector_A1   = zeros(SSdim_ID-(it-1), 2)
        vector_A2   = zeros(SSdim_ID-(it-1), 1)

        vector_b1   = zeros(SSdim_ID-(it-1), 1)
        vector_b2   = zeros(SSdim_ID-(it-1), 1)

        Counter_ID = 1
        for ii = 1:(it-1)
            for jj = 1:(time[ii]-1)
                vector_A1[Counter_ID,:] = [ SS[1,jj,ii], SS[2,jj,ii] ]
                vector_b1[Counter_ID,:] = [ SS[1,jj+1,ii] - B[1] * OldU[1,jj,ii]]
                
                vector_A2[Counter_ID,:] = [ SS[2,jj,ii] ]
                vector_b2[Counter_ID,:] = [ SS[2,jj+1,ii] - B[2] * OldU[1,jj,ii]]

                Counter_ID = Counter_ID + 1
            end
        end
        for ii = 1:t
            vector_A1[Counter_ID,:] = [ x_LMPC[1,ii], x_LMPC[2,ii] ]
            vector_b1[Counter_ID,:] = [ x_LMPC[1,ii+1] - B[1] * u_LMPC[1,ii]]
            
            vector_A2[Counter_ID,:] = [ x_LMPC[2,ii] ]
            vector_b2[Counter_ID,:] = [ x_LMPC[2,ii+1] - B[2] * u_LMPC[1,ii]]

            Counter_ID = Counter_ID + 1
        end


        Matrix1 = vector_A1'*vector_A1
        Matrix2 = vector_A2'*vector_A2


        Row1    = inv(Matrix1) * vector_A1' * vector_b1
        Row2    = inv(Matrix2) * vector_A2' * vector_b2

        MeanEstimate = zeros(2,2)
        MeanEstimate[1,:] = Row1
        MeanEstimate[2,2] = Row2[1]

        MSE1 = 0;
        MSE2 = 0;
        for i = 1:(SSdim_ID-it)
            MSE1 = MSE1 + ( vector_b1[i,:] - *(vector_A1[i,:], Row1) )^2
            MSE2 = MSE2 + ( vector_b2[i,:] - vector_A2[i,1]* Row2[1] )^2
        end
        MSE1 = MSE1^(0.5)/(SSdim_ID-it-1)
        MSE2 = MSE2^(0.5)/(SSdim_ID-it-1)

        MSE = zeros(2)
        MSE[1] = MSE1[1]
        MSE[2] = MSE2[1]

        println("Estimte is ", MeanEstimate)
        t=t+1
    end
    # ========================================================================================================
    # ==================================== Back to the iterations loop =======================================
    # ========================================================================================================

    # Now post process the data after the LMPC has converged
    time[it] = t

    # Add data to SS and Q function
    OldU[:, 1:time[it]-1, it] = u_LMPC[:,1:time[it]-1]
    SS[:, 1:time[it], it]   = x_LMPC[:,1:time[it]]
    Qfun[:, 1:time[it], it] = ComputeCost(x_LMPC[:,1:time[it]], u_LMPC[:,1:time[it]], LMPCparams)
    OptQ[:, 1:time[it], it] = ComputeCost(OptX[:,1:time[it],it], OldU[:,1:time[it],it], LMPCparams)



    Difference = Qfun[1,1,it-1]-Qfun[1,1,it]


    SSdim_ID = sum(time)
    vector_A1   = zeros(SSdim_ID-it, 2)
    vector_A2   = zeros(SSdim_ID-it, 1)

    vector_b1   = zeros(SSdim_ID-it, 1)
    vector_b2   = zeros(SSdim_ID-it, 1)

    Counter_ID = 1
    for ii = 1:it
        for jj = 1:(time[ii]-1)
            vector_A1[Counter_ID,:] = [ SS[1,jj,ii], SS[2,jj,ii] ]
            vector_b1[Counter_ID,:] = [ SS[1,jj+1,ii] - B[1] * OldU[1,jj,ii]]
            
            vector_A2[Counter_ID,:] = [ SS[2,jj,ii] ]
            vector_b2[Counter_ID,:] = [ SS[2,jj+1,ii] - B[2] * OldU[1,jj,ii]]

            Counter_ID = Counter_ID + 1
        end
    end

    Matrix1 = vector_A1'*vector_A1
    Matrix2 = vector_A2'*vector_A2


    Row1    = inv(Matrix1) * vector_A1' * vector_b1
    Row2    = inv(Matrix2) * vector_A2' * vector_b2

    MeanEstimate = zeros(2,2)
    MeanEstimate[1,:] = Row1
    MeanEstimate[2,2] = Row2[1]

    MSE1 = 0;
    MSE2 = 0;
    for i = 1:(SSdim_ID-it)
        MSE1 = MSE1 + ( vector_b1[i,:] - *(vector_A1[i,:], Row1) )^2
        MSE2 = MSE2 + ( vector_b2[i,:] - vector_A2[i,1]* Row2[1] )^2
    end
    MSE1 = MSE1^(0.5)/(SSdim_ID-it-1)
    MSE2 = MSE2^(0.5)/(SSdim_ID-it-1)

    MSE = zeros(2)
    MSE[1] = MSE1[1]
    MSE[2] = MSE2[1]

    it = it + 1

end
it = it - 1

figure()
hold(1)
for i = 1:it-1
    plot(SS[1,:,i]',SS[2,:,i]', "-ks")
end
plot(x_LMPC[1,:]',x_LMPC[2,:]', "-ro")
plot(OptX[1,:,it]',OptX[2,:,it]', "-g*")


figure()
index = it
plot(SS[1,:,index]',SS[2,:,index]', "-ko")
plot(OptX[1,:,index]',OptX[2,:,index]', "-g*")


figure()
index = it-1
plot(SS[1,:,index]',SS[2,:,index]', "-ko")
plot(OptX[1,:,index]',OptX[2,:,index]', "-g*")

for i = 1:it-1
    println(i,"-th itearion cost LMPC; ", Qfun[1,1,i])
    println(i,"-th itearion cost; ", OptQ[1,1,i])
end
# for i = 2:it
#     figure()
#     plot(x_feasible[1,:]',x_feasible[2,:]', "-ro")
#     hold(1)

#     j = 1
#     plot(SS[1, 1:time[j], j]', SS[2, 1:time[j], j]', "-go", label="Safe Set")
#     for j = 2:(i-1)
#         plot(SS[1, 1:time[j], j]', SS[2, 1:time[j], j]', "-go")
#     end

#     plot(SS[1, 1:time[i], i]', SS[2, 1:time[i], i]', "-ko", label="Trajectory Reference System")
#     plot(x_real[1, 1:time[i], i]', x_real[2, 1:time[i], i]', "-r*", label="Trajectory Real System" )
    
#     grid(1)
#     LMPCIteration = i -1
#     title("Closed-loop trajectory at iteration $LMPCIteration")
#     axis("equal")
#     xlabel(L"$x_1$", size=24)
#     ylabel(L"$x_2$", size=24)
#     legend(loc="lower right",fancybox="true")
# end

# figure()
# hold(1)
# i = 2
# vec1 = SS[3, 1:time[i], i].^2 + SS[4, 1:time[i], i].^2
# println("Here ",vec)
# plot(1:time[i], vec1[:] , "-ro", label="1st Iteration")

# i = it
# vecss = SS[3, 1:time[i], i].^2 + SS[4, 1:time[i], i].^2
# println("Here ",vec)
# plot(1:time[i], vecss[:] , "-g*", label="Steady State")

# grid(1)
# title("LMPC Steady State")
# axis("equal")
# xlabel(L"Time step", size=24)
# ylabel(L"$||x_2||_2^2$", size=24)
# legend()
