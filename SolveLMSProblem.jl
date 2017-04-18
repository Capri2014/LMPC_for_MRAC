function solveLMSProblem(mdl::LMS_Model,x_feas::Array{Float64,2})

    # Load Parameters
    sol_status::Symbol

    # Update current initial condition, curvature and previous input
    setvalue(mdl.x_feas,x_feas)

    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)

    #LMPCSol.x    = getvalue(mdl.k_til)

    println("Solved, status = $sol_status")
    return getvalue(mdl.k_til), getvalue(mdl.Dummy)
end