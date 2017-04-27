type NominalLMPC_Model
    mdl::JuMP.Model

    x0::Array{JuMP.NonlinearParameter,1}
    SS::Array{JuMP.NonlinearParameter,2}
    Qfun::Array{JuMP.NonlinearParameter,1}
    Mean::Array{JuMP.NonlinearParameter,2}

    x_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}
    lamb::Array{JuMP.Variable,2}

    state_cost::JuMP.NonlinearExpression
    input_cost::JuMP.NonlinearExpression
    termi_cost::JuMP.NonlinearExpression

    function NominalLMPC_Model(LMPCparams::TypeLMPCparams,SystemParams::TypeSystemParams, SSdim::Int64, Mean::Array{Float64,2})
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
        @variable( mdl, lamb[1:SSdim,1])


        @NLparameter(mdl, x0[1:n] == 0)
        @NLparameter(mdl, SS[1:n,1:SSdim] == 0)
        @NLparameter(mdl, Qfun[1:SSdim] == 0)
        @NLparameter(mdl, Mean[1:2,1:2] == 0)

        # System dynamics
        @NLconstraint(mdl, [i=1:n], x_Ol[i,1] == x0[i])         # initial condition
        println("Initializing model...")

        # System dynamics

        for i=1:N           
            @NLconstraint(mdl, x_Ol[1,i+1] == Mean[1,1] * x_Ol[1,i] + Mean[1,2] * x_Ol[2,i])
            @NLconstraint(mdl, x_Ol[2,i+1] ==    0                  + Mean[2,2] * x_Ol[2,i] + B[2]*u_Ol[1,i])
        end

        # Constratints Related with the LMPC
        for j=1:SSdim
            setlowerbound(lamb[j,1], 0)
        end
        @NLconstraint(mdl,sum{lamb[j,1], j = 1:SSdim} == 1)


        @NLconstraint(mdl, x_Ol[1,N+1] == sum{SS[1,k] * lamb[k,1], k = 1:SSdim})        
        @NLconstraint(mdl, x_Ol[2,N+1] == sum{SS[2,k] * lamb[k,1], k = 1:SSdim})         
       
        # Cost definitions
        # State cost
        # ---------------------------------
        @NLexpression(mdl, state_cost, sum{sum{ (Q[j,j]  * x_Ol[j,i])^2  , i=1:N},j=1:2})

        # Control Input cost
        @NLexpression(mdl, input_cost, sum{ (R[1,1] * u_Ol[1,i])^2, i=1:N})

        # Control Input cost
        @NLexpression(mdl, termi_cost, sum{ Qfun[j] * lamb[j,1] ,j=1:SSdim})


        # Objective function
        @NLobjective(mdl, Min, state_cost + input_cost + termi_cost)

        # First solve
        #sol_stat=solve(mdl)
        #println("Finished solve 1: $sol_stat")
        #sol_stat=solve(mdl)
        #println("Finished solve 2: $sol_stat")
        
        m.mdl  = mdl
        m.x0   = x0
        m.Qfun = Qfun
        m.SS   = SS
        m.Mean   = Mean

        m.x_Ol = x_Ol
        m.u_Ol = u_Ol

        m.state_cost = state_cost
        m.input_cost = input_cost
        m.termi_cost = termi_cost
        
        return m
    end
end