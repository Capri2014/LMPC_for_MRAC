function Feasible_Traj(SystemParams::TypeSystemParams, x0::Array{Float64,1}, K_r::Array{Float64,1})

    Ar = SystemParams.Ar
    B  = SystemParams.B
    n  = SystemParams.n
    d  = SystemParams.d

    Points = 100
    PointsSysID = 80
    
    x_real = zeros(n, Points+1)
    x_real[:,1] = x0
    u_real = zeros(d, Points)

    
    x_feasible = zeros(2, Points+1)
    u_feasible = zeros(1, Points)

    vector_A1   = zeros(PointsSysID, 2)
    vector_A2   = zeros(PointsSysID, 1)

    vector_b1   = zeros(PointsSysID, 1)
    vector_b2   = zeros(PointsSysID, 1)

    for i = 1:Points
        u_real[:,i]   = dot(K_r, x_real[:,i])
        noise = 1*[2*randn(); 3*randn()]
        x_real[:,i+1] = Ar * x_real[:,i] + B * u_real[1,i] + noise
        u_apply = B * u_real[1,i]

        x_feasible[:,i] = x_real[:,i]
        u_feasible[:,i] = u_real[1,i]
        if (i > 1) && (i<PointsSysID+1)
            vector_A1[i-1,:] = [ x_feasible[1,i-1], x_feasible[2,i-1] ]
            vector_b1[i-1,:] = [ x_feasible[1,i] - B[1] * u_feasible[1,i-1]]
            
            vector_A2[i-1,:] = [ x_feasible[2,i-1]]
            vector_b2[i-1,:] = [ x_feasible[2,i] - B[2] * u_feasible[1,i-1]]
        end

    end
    i = Points + 1
    x_feasible[:,i] = x_real[:,i]

Matrix1 = vector_A1'*vector_A1
Matrix2 = vector_A2'*vector_A2


Row1    = inv(Matrix1) * vector_A1' * vector_b1
Row2    = inv(Matrix2) * vector_A2' * vector_b2

MeanEstimate = zeros(2,2)
MeanEstimate[1,:] = Row1
MeanEstimate[2,2] = Row2[1]

MSE1 = 0;
MSE2 = 0;
for i = 1:PointsSysID
    
    MSE1 = MSE1 + ( vector_b1[i,:] - *(vector_A1[i,:], Row1) )^2
    MSE2 = MSE2 + ( vector_b2[i,:] - vector_A2[i,1]* Row2[1] )^2
end
MSE1 = MSE1^(0.5)/(PointsSysID-1)
MSE2 = MSE2^(0.5)/(PointsSysID-1)

MSE = zeros(2)
MSE[1] = MSE1[1]
MSE[2] = MSE2[1]

return x_feasible, u_feasible, MeanEstimate, MSE


end