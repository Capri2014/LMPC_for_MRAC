function ComputeCost(x::Array{Float64,2}, u::Array{Float64,2}, LMPCparams::TypeLMPCparams, SystemParams::TypeSystemParams)

    Points = size(x)[2] 
    Qfun = zeros(Points)
    xF  = SystemParams.xF
    dt  = SystemParams.dt
    rho = SystemParams.rho
    g   = SystemParams.g

    for i = 1:Points
        index = Points - i + 1
        if i == 1
            v = x[1, index]
            s = x[2, index]
            ov= x[3, index]
            os= x[4, index]

            otheta = 3*sin( (os - xF[2]) / xF[2] * 2 * 3.14 )

	    Norm2 = v^2 + ( s - xF[2] )^2 + ov^2 +( os - xF[4] )^2 
	    Qfun[index] = 10*Norm2 / ((10*Norm2)^2 + 1)^0.5 +  ( (v-ov)/dt + rho*ov^2 + g*sin(otheta))^2 
        else
            v = x[1, index]
            s = x[2, index]
            ov= x[3, index]
            os= x[4, index]

            otheta = 3*sin( (os - xF[2]) / xF[2] * 2 * 3.14 )

            Norm2 = v^2 + ( s - xF[2] )^2 + ov^2 +( os - xF[4] )^2 
	    Qfun[index] = Qfun[index+1] + 10*Norm2/( (10*Norm2)^2  + 1)^0.5 + ( (v-ov)/dt + rho*ov^2 + g*sin(otheta))^2  
        end
    end

    return Qfun
    
end
