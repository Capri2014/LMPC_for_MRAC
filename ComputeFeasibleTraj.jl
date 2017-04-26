function Feasible_Traj(SystemParams::TypeSystemParams, x0::Array{Float64,1}, K_r::Array{Float64,2})
include("LMS_Model.jl")
include("SolveLMSProblem.jl")


    Ad = SystemParams.Ad
    Ar = SystemParams.Ar
    B  = SystemParams.B
    n  = SystemParams.n
    d  = SystemParams.d

    Points = 180
    PointsSysID = 5
    
    x_real = zeros(n, Points+1)
    x_real[:,1] = x0
    u_real = zeros(d, Points)


    x_desi = zeros(n, Points+1)
    x_desi[:,1] = x0


    K_tilda = K_r - (Ad - Ar)

    
    x_feasible = zeros(7, Points+1)
    u_feasible = zeros(5, Points)

    vector_b   = zeros(PointsSysID, 1)
    vector_A   = zeros(PointsSysID, 1)

    vector_b1   = zeros(PointsSysID, 1)
    vector_A1   = zeros(PointsSysID, 2)

    for i = 1:Points
        u_real[:,i]   = K_r * x_real[:,i]
        x_real[:,i+1] = Ar * x_real[:,i] + u_real[:,i]


        x_desi[:,i+1] = Ad * x_desi[:,i]

        x_feasible[1:2,i] = x_real[:,i]
        x_feasible[3:4,i] = x_real[:,i] - x_desi[:,i]
        x_feasible[5,i]   = K_tilda[1,1]
        x_feasible[6,i]   = K_tilda[1,2]
        x_feasible[7,i]   = K_tilda[2,2]

        if (i > 1) && (i<PointsSysID+1)
            vector_A1[i-1,:] = [x_feasible[1,i-1] + x_feasible[3,i-1], x_feasible[2,i-1] + x_feasible[4,i-1]]
            vector_b1[i-1,:] = x_feasible[3,i] - Ad[1,1] * x_feasible[3,i-1]-Ad[1,2] * x_feasible[4,i-1]
            
            vector_A[i-1,1] = x_feasible[2,i-1] + x_feasible[4,i-1]
            vector_b[i-1,1] = x_feasible[4,i] - Ad[2,1] * x_feasible[3,i-1]-Ad[2,2] * x_feasible[4,i-1]

        end

        u_feasible[1:2,i] = u_real[:,i]
        u_feasible[3,i] = 0
        u_feasible[4,i] = 0
        u_feasible[5,i] = 0

    end
    i = Points + 1
    
    x_feasible[1:2,i] = x_real[:,i]
    x_feasible[3:4,i] = x_real[:,i] - x_desi[:,i]
    x_feasible[5,i]   = K_tilda[1,1]
    x_feasible[6,i]   = K_tilda[1,2]
    x_feasible[7,i]   = K_tilda[2,2]

Matrix = vector_A1'*vector_A1
K12 = inv(Matrix) * vector_A1' * vector_b1

Matrix = vector_A'*vector_A
K22 = inv(Matrix) * vector_A' * vector_b

for i = 1:Points
    x_feasible[5,i]   = K12[1]
    x_feasible[6,i]   = K12[2]
    x_feasible[7,i]   = K22[1]

end
i = Points + 1

x_feasible[5,i]   = K12[1]
x_feasible[6,i]   = K12[2]
x_feasible[7,i]   = K22[1]

MatrixPlot = [K12[1] K12[2]; 0 K22]
println("Estimated K_tilda", MatrixPlot)
println("Real K_tilda", K_tilda)

return x_feasible, u_feasible


end