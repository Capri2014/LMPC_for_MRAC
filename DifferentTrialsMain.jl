using JuMP
using Ipopt
using JLD
using PyPlot
using JLD
using PyPlot

include("DifferentTrials.jl")


Val = DifferentTrials()

BetterOpr = 0
WorseOpr = 0

for i = 1:2000
    close("all")
    println(i)
    Dummy = DifferentTrials()
    if Dummy < 0
        BetterOpr = BetterOpr + 1
    else
        WorseOpr  = WorseOpr + 1
    end
    Val = Val + Dummy
    println(Val, " BetterOpr: ", BetterOpr, "  WorseOpr: ", WorseOpr)
end

Mean = Val/(50+1)

println(Val, " BetterOpr: ", BetterOpr, "  WorseOpr: ", WorseOpr)
