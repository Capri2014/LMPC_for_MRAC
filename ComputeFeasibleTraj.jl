function Feasible_Traj(SystemParams::TypeSystemParams, x0::Array{Float64,1}, K_r::Array{Float64,2})
include("LMS_Model.jl")
include("SolveLMSProblem.jl")


    Ad = SystemParams.Ad
    Ar = SystemParams.Ar
    B  = SystemParams.B
    n  = SystemParams.n
    d  = SystemParams.d

    Points = 150
    x_real = zeros(n, Points+1)
    x_real[:,1] = x0
    u_real = zeros(d, Points)


    x_desi = zeros(n, Points+1)
    x_desi[:,1] = x0


    K_tilda = K_r - (Ad - Ar)

    
    x_feasible = zeros(7, Points+1)
    u_feasible = zeros(5, Points)

    for i = 1:Points
        u_real[:,i]   = K_r * x_real[:,i]
        x_real[:,i+1] = Ar * x_real[:,i] + B * u_real[:,i]


        x_desi[:,i+1] = Ad * x_desi[:,i]

        x_feasible[1:2,i] = x_real[:,i]
        x_feasible[3:4,i] = x_real[:,i] - x_desi[:,i]
        x_feasible[5,i]   = K_tilda[1,1]
        x_feasible[6,i]   = K_tilda[1,2]
        x_feasible[7,i]   = K_tilda[2,2]

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



mdl_LMS     = LMS_Model(SystemParams, Points+1)
println("Here")

k_tilda_LMS, Dummy = solveLMSProblem(mdl_LMS, x_feasible[1:4,:])

println("Dummy is ", Dummy)

println("K_tilda is ", K_tilda)

println("K_tilda optimization is ", k_tilda_LMS)

for i = 1:Points
    x_feasible[5,i]   = k_tilda_LMS[1]
    x_feasible[6,i]   = k_tilda_LMS[2]
    x_feasible[7,i]   = k_tilda_LMS[3]

end
i = Points + 1

x_feasible[5,i]   = k_tilda_LMS[1]
x_feasible[6,i]   = k_tilda_LMS[2]
x_feasible[7,i]   = k_tilda_LMS[3]

return x_feasible, u_feasible


end