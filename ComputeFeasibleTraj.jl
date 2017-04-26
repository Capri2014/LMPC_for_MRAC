function Feasible_Traj(SystemParams::TypeSystemParams, x0::Array{Float64,1}, K_r::Array{Float64,1})

    Ar = SystemParams.Ar
    B  = SystemParams.B
    n  = SystemParams.n
    d  = SystemParams.d

    Points = 180
    PointsSysID = 150
    
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
        noise =0.0*[2*randn(); 3*randn()]
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

MSE2 = 0;
for i = 1:PointsSysID
    MSE2 = MSE2 + ( vector_b2[i,:] - vector_A2[i,1]* Row2[1] )^2
    println(vector_b2[i,:] - vector_A2[i,1]* Row2[1])
end
MSE2 = MSE2/(PointsSysID-1)

println(Row1)
println(Row2)
println("MES",MSE2)

# for i = 1:Points
#     x_feasible[5,i]   = K12[1]
#     x_feasible[6,i]   = K12[2]
#     x_feasible[7,i]   = K22[1]

# end
# i = Points + 1

# x_feasible[5,i]   = K12[1]
# x_feasible[6,i]   = K12[2]
# x_feasible[7,i]   = K22[1]

# MatrixPlot = [K12[1] K12[2]; 0 K22]
# println("Estimated K_tilda", MatrixPlot)
# println("Real K_tilda", K_tilda)

return x_feasible, u_feasible


end