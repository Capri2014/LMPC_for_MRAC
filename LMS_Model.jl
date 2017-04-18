type LMS_Model
    mdl::JuMP.Model

    x_feas::Array{JuMP.NonlinearParameter,2}

    Dummy::Array{JuMP.Variable,2}
    k_til::Array{JuMP.Variable,1}

    cost_l2::JuMP.NonlinearExpression

    function LMS_Model(SystemParams::TypeSystemParams, Npoints::Int64)
        println("Starting creation ==================================<=================")
        m = new()
        A          = SystemParams.Ad
        B          = SystemParams.B

        # Create Model
        mdl = Model(solver = IpoptSolver(print_level=1))

        # Create variables (these are going to be optimized)
        @variable( mdl, Dummy[1:2,1:Npoints-1]) 
        @variable( mdl, k_til[1:3]) 


        @NLparameter(mdl, x_feas[1:4,1:Npoints] == 0)


        # System dynamics
        println("Initializing model...")

        # System dynamics

        for i=1:(Npoints-1)
            @NLconstraint(mdl, Dummy[1,i] == - x_feas[3,i+1] + sum{A[1,k] * x_feas[k+2,i], k=1:2}
                         + k_til[1]*(x_feas[1,i] + x_feas[3,i]) + k_til[2]*(x_feas[2,i] + x_feas[4,i]) )
           
            @NLconstraint(mdl, Dummy[2,i] == - x_feas[4,i+1] + sum{A[2,k] * x_feas[k+2,i], k=1:2} 
                         + k_til[3]*(x_feas[2,i] + x_feas[4,i]) )

        end
         # Cost definitions
        # State cost
        # ---------------------------------
        @NLexpression(mdl, cost_l2, sum{sum{ (Dummy[j,i])^2, i=1:(Npoints-1)},j=1:2})

        # Objective function
        @NLobjective(mdl, Min, cost_l2)

        # First solve
        #sol_stat=solve(mdl)
        #println("Finished solve 1: $sol_stat")
        #sol_stat=solve(mdl)
        #println("Finished solve 2: $sol_stat")
        
        m.mdl     = mdl
        m.x_feas  = x_feas
        m.k_til   = k_til
        m.Dummy   = Dummy
        m.cost_l2 = cost_l2
        
        return m
    end
end