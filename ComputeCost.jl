function ComputeCost(x::Array{Float64,2}, u::Array{Float64,2}, LMPCparams::TypeLMPCparams)

    N           = LMPCparams.N
    Q           = LMPCparams.Q
    Qe          = LMPCparams.Qe
    R           = LMPCparams.R

    Points = size(x)[2] 
    Qfun = zeros(Points)

    for i = 1:Points

        index = Points -i+1
        if i == 1
            Qfun[index] = Q[1,1] * x[1, index]^2 + Q[2,2] * x[2, index]^2 + Qe[1,1] * x[3, index]^2 + Qe[2,2] * x[4, index]^2  
        else

            Qfun[index] = Qfun[index+1] + Q[1,1]  * x[1, index]^2 + Q[2,2]  * x[2, index]^2 +
                                          Qe[1,1] * x[3, index]^2 + Qe[2,2] * x[4, index]^2 +
                                          R[1,1]  * u[1, index]^2 + R[2,2]  * u[2, index]^2
        end
    end

    return Qfun
    
end