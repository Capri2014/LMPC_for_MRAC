using JuMP
using Ipopt
using JLD
using PyPlot
using JLD
using PyPlot

include("DifferentTrials.jl")


Val = DifferentTrials()

for i = 1:200
    close("all")
    println(i)
    Dummy = DifferentTrials()
    Val = Val + Dummy
    println(Val)
end

Mean = Val/(50+1)

println(Val)
