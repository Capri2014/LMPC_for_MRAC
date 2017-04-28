function solveLMPCProblem(mdl::LMPC_Model,LMPCSol::TypeLMPCSol,xCurr::Array{Float64,1},Mean::Array{Float64,2},Variance::Array{Float64,1})

    # Load Parameters
    sol_status::Symbol

    # Update current initial condition, curvature and previous input
    setvalue(mdl.x0,xCurr)
    setvalue(mdl.Mean,Mean)
    setvalue(mdl.Variance,Variance)

    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)

    LMPCSol.x    = getvalue(mdl.x_Ol)
    LMPCSol.u    = getvalue(mdl.u_Ol)
    LMPCSol.a    = getvalue(mdl.a_Ol)

    LMPCSol.cost = getvalue(mdl.state_cost) + getvalue(mdl.input_cost)

    # println("LMPC Solved, status = $sol_status")
    nothing
end