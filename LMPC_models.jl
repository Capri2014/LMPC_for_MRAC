type LMPC_Model
    mdl::JuMP.Model

    x0::Array{JuMP.NonlinearParameter,1}
    Mean::Array{JuMP.NonlinearParameter,2}
    V::Array{JuMP.NonlinearParameter,2}
    Variance::Array{JuMP.NonlinearParameter,1}
    beta::Array{JuMP.NonlinearParameter,1}

    x_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}
    a_Ol::Array{JuMP.Variable,1}
    lamb::Array{JuMP.Variable,2}

    state_cost::JuMP.NonlinearExpression
    input_cost::JuMP.NonlinearExpression

    function LMPC_Model(LMPCparams::TypeLMPCparams,SystemParams::TypeSystemParams)
        println("Starting creation")
        m = new()
        B          = SystemParams.B
        n          = 2
        d          = 1

        N           = LMPCparams.N
        Q           = LMPCparams.Q
        R           = LMPCparams.R

        # Create Model
        mdl = Model(solver = IpoptSolver(print_level=0))

        # Create variables (these are going to be optimized)
        @variable( mdl, x_Ol[1:n,1:(N+1)]) 
        @variable( mdl, u_Ol[1:d,1:N])
        @variable( mdl, a_Ol[1:6])
        
        # setvalue(a_Ol[1],0)
        # setvalue(a_Ol[2],0)
        # setvalue(a_Ol[3],0)
        # setvalue(a_Ol[4],0)
        # setvalue(a_Ol[5],0)
        # setvalue(a_Ol[6],0)

        @NLparameter(mdl, x0[1:n] == 0)
        @NLparameter(mdl, Mean[1:2,1:3] == 0)
        @NLparameter(mdl, V[1:3,1:3] == 0)
        @NLparameter(mdl, beta[1:2] == 0)
        @NLparameter(mdl, Variance[1:2] == 0)





        # System dynamics
        @NLconstraint(mdl, [i=1:n], x_Ol[i,1] == x0[i])         # initial condition
        println("Initializing model...")

        # System dynamics

        for i=1:N           
            @NLconstraint(mdl, x_Ol[1,i+1] == a_Ol[1] * x_Ol[1,i] + a_Ol[2] * x_Ol[2,i] + a_Ol[5] * u_Ol[1,i])
            @NLconstraint(mdl, x_Ol[2,i+1] == a_Ol[3] * x_Ol[1,i] + a_Ol[4] * x_Ol[2,i] + a_Ol[6] * u_Ol[1,i])
        end

        @NLconstraint(mdl, a_Ol[1]^2 + a_Ol[2]^2 + a_Ol[3]^2 + 
                           a_Ol[4]^2 + a_Ol[5]^2 + a_Ol[6]^2 <= 100)

        # @NLconstraint(mdl, V[1,1] *( (a_Ol[1]-Mean[1,1])^2 + (a_Ol[3]-Mean[2,1])^2 )   
        #                   +V[2,2] *( (a_Ol[2]-Mean[1,2])^2 + (a_Ol[4]-Mean[2,2])^2 )   
        #                   +V[3,3] *( (a_Ol[5]-Mean[1,3])^2 + (a_Ol[6]-Mean[2,3])^2 )
        #                 +2*V[1,2] *( (a_Ol[1]-Mean[1,1]) * (a_Ol[2]-Mean[1,2]) + (a_Ol[4]-Mean[2,2]) * (a_Ol[3]-Mean[2,1]) ) 
        #                 +2*V[1,3] *( (a_Ol[1]-Mean[1,1]) * (a_Ol[5]-Mean[1,3]) + (a_Ol[3]-Mean[2,1]) * (a_Ol[6]-Mean[2,3]) )
        #                 +2*V[2,3] *( (a_Ol[3]-Mean[1,2]) * (a_Ol[5]-Mean[1,3]) + (a_Ol[4]-Mean[2,2]) * (a_Ol[6]-Mean[2,3]) ) <= beta[1])
        
       
        for i=1:N            
            @NLconstraint(mdl, (a_Ol[1] - Mean[1,1]) * x_Ol[1,i] + (a_Ol[2] - Mean[1,2]) * x_Ol[2,i] + (a_Ol[5] - Mean[1,3]) * u_Ol[1,i] >= -0.1*Variance[1])
            @NLconstraint(mdl, (a_Ol[1] - Mean[1,1]) * x_Ol[1,i] + (a_Ol[2] - Mean[1,2]) * x_Ol[2,i] + (a_Ol[5] - Mean[1,3]) * u_Ol[1,i] <=  0.1*Variance[1])
            @NLconstraint(mdl, (a_Ol[3] - Mean[2,1]) * x_Ol[1,i] + (a_Ol[4] - Mean[2,2]) * x_Ol[2,i] + (a_Ol[6] - Mean[2,3]) * u_Ol[1,i] >= -0.1*Variance[2])
            @NLconstraint(mdl, (a_Ol[3] - Mean[2,1]) * x_Ol[1,i] + (a_Ol[4] - Mean[2,2]) * x_Ol[2,i] + (a_Ol[6] - Mean[2,3]) * u_Ol[1,i] <=  0.1*Variance[2])
        end

        # Constratints Related with the LMPC


        # @NLconstraint(mdl, x_Ol[1,N+1] == sum{SS[1,k] * lamb[k,1], k = 1:SSdim})        
        # @NLconstraint(mdl, x_Ol[2,N+1] == sum{SS[2,k] * lamb[k,1], k = 1:SSdim})         
       
        # Cost definitions
        # State cost
        # ---------------------------------
        @NLexpression(mdl, state_cost, sum{sum{ (Q[j,j]  * x_Ol[j,i])^2  , i=1:N},j=1:2})

        # Control Input cost
        @NLexpression(mdl, input_cost, sum{ (R[1,1] * u_Ol[1,i])^2, i=1:N})
                                            # + 0.00001*(a_Ol[1] - Mean[1,1])^2
                                            # + 0.00001*(a_Ol[2] - Mean[1,2])^2
                                            # + 0.00001*(a_Ol[3] - Mean[2,1])^2
                                            # + 0.00001*(a_Ol[4] - Mean[2,2])^2
                                            # + 0.00001*(a_Ol[5] - Mean[1,3])^2
                                            # + 0.00001*(a_Ol[6] - Mean[2,3])^2)

        # Control Input cost
        # @NLexpression(mdl, termi_cost, sum{ Qfun[j] * lamb[j,1] ,j=1:SSdim})


        # Objective function
        @NLobjective(mdl, Min, state_cost + input_cost )

        # First solve
        #sol_stat=solve(mdl)
        #println("Finished solve 1: $sol_stat")
        #sol_stat=solve(mdl)
        #println("Finished solve 2: $sol_stat")
        
        m.mdl  = mdl
        m.x0   = x0
        m.Variance = Variance
        m.Mean   = Mean
        m.beta = beta

        m.x_Ol = x_Ol
        m.V    = V

        m.u_Ol = u_Ol
        m.a_Ol = a_Ol

        m.state_cost = state_cost
        m.input_cost = input_cost
        
        return m
    end
end