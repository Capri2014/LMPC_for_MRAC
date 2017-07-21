function Feasible_Traj(SystemParams::TypeSystemParams, x0::Array{Float64,1})
include("LMS_Model.jl")
include("SolveLMSProblem.jl")


    g   = SystemParams.g
    dt  = SystemParams.dt
    xF  = SystemParams.xF
    rho = SystemParams.rho


    Points = 200
    
    x_feasible = zeros(4, Points+1)
    u_feasible = zeros(1, Points)


    x_feasible[:,1] = x0
    u_feasible[:,1] = 5

    for i = 1:Points
	    teta = 3*sin( (x_feasible[2,i] - xF[2]) / xF[2] * 2 * 3.14 )
	
	# Logic to compute input
	if x_feasible[1,i] > 10
		u_feasible[1,i] = 0
	else
		u_feasible[1,i] = u_feasible[1,1]  + rho * x_feasible[1,i] + g * sin(teta)
	end
	
	Next_s = x_feasible[2,i] + dt * x_feasible[1,i]
	if (xF[2] - Next_s ) < 5
		u_feasible[1,i] = -0.7* u_feasible[1,1]  + rho * x_feasible[1,i] + g * sin(teta)
	end
	if (xF[2] - Next_s ) < 0.5
		v_desired       = (xF[2] - Next_s)/dt
		u_feasible[1,i] = ( v_desired - x_feasible[1,i]  )/dt + rho*x_feasible[1,i]^ 2 + g*sin(teta) 
	end

	if xF[2] == Next_s
		u_feasible[1,i] = ( 0 - x_feasible[1,i]  )/dt + rho* x_feasible[1,i]^2 + g*sin(teta) 
	end

	# Propagate input forward
	x_feasible[1,i+1] = x_feasible[1,i] + dt*( u_feasible[1,i]  - rho* x_feasible[1,i]^2 - g*sin( teta ) )
	x_feasible[2,i+1] = x_feasible[2,i] + dt*( x_feasible[1,i] )
	x_feasible[3,i+1] = x_feasible[1,i]
	x_feasible[4,i+1] = x_feasible[2,i]


    end
    
    return x_feasible, u_feasible

end
