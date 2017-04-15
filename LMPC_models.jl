type LMPC_Model
    mdl::JuMP.Model

    x0::Array{JuMP.NonlinearParameter,1}

    x_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}

    state_cost::JuMP.NonlinearExpression
    input_cost::JuMP.NonlinearExpression

    function LMPC_Model(LMPCparams::TypeLMPCparams,SystemParams::TypeSystemParams)
        println("Starting creation")
        m = new()
        A          = SystemParams.A
        B          = SystemParams.B
        n          = SystemParams.n
        d          = SystemParams.d

        N           = LMPCparams.N
        Q           = LMPCparams.Q
        R           = LMPCparams.R

        # Create Model
        mdl = Model(solver = IpoptSolver(print_level=0))

        # Create variables (these are going to be optimized)
        @variable( mdl, x_Ol[1:(N+1),1:n]) 
        @variable( mdl, u_Ol[1:N,1:d])

        @NLparameter(mdl, x0[i=1:n] == 0)


        # System dynamics
        @NLconstraint(mdl, [i=1:n], x_Ol[1,i] == x0[i])         # initial condition
        println("Initializing model...")

        # System dynamics
        for i=1:N
            for j=1:2
                    @NLconstraint(mdl, x_Ol[i+1,j] == sum{A[j,k] * x_Ol[i,k] + B[j,k] * u_Ol[i,k], k=1:2})
            end
        end

        # Cost definitions
        # State cost
        # ---------------------------------
        @NLexpression(mdl, state_cost, sum{sum{ (Q[j,j] * x_Ol[i,j])^2, i=1:N},j=1:n})

        # Control Input cost
        @NLexpression(mdl, input_cost, sum{sum{ (R[j,j] * u_Ol[i,j])^2, i=1:N},j=1:n})


        # Objective function
        @NLobjective(mdl, Min, state_cost + input_cost)

        # create first artificial solution (for warm start)
        # for i=1:N+1
        #     setvalue(z_Ol[i,:],[(i-1)*dt*v_ref 0 0 v_ref 0])
        # end
        # for i=1:N
        #     setvalue(u_Ol[i,:],[0.15 0])
        # end

        # First solve
        sol_stat=solve(mdl)
        println("Finished solve 1: $sol_stat")
        sol_stat=solve(mdl)
        println("Finished solve 2: $sol_stat")
        
        m.mdl = mdl
        m.x0 = x0
        m.x_Ol = x_Ol
        m.u_Ol = u_Ol
        m.state_cost = state_cost
        m.input_cost = input_cost
        return m
    end
end