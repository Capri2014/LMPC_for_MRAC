function solveLMPCProblem(mdl::LMPC_Model,LMPCSol::TypeLMPCSol,xCurr::Array{Float64,1},ConvSS::Array{Float64,2},ConvQfun::Array{Float64,1})

    # Load Parameters
    sol_status::Symbol

    # Update current initial condition, curvature and previous input
    setvalue(mdl.x0,xCurr)
    setvalue(mdl.SS,ConvSS)
    setvalue(mdl.Qfun,ConvQfun)

    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)

    LMPCSol.x    = getvalue(mdl.x_Ol)
    LMPCSol.u    = getvalue(mdl.u_Ol)

    LMPCSol.cost = getvalue(mdl.state_cost) + getvalue(mdl.termi_cost)

    println("Solved, status = $sol_status")
    nothing
end
