using JuMP
using Ipopt
using JLD
using PyPlot
using JLD
using PyPlot

include("DifferentTrials.jl")


Val = DifferentTrials()

for i = 1:2000
    println(i)
    Dummy = DifferentTrials()
    Val = Val + Dummy
    println("Cumulate", Val)
end

Mean = Val/(50+1)

println("Mean",Mean)
