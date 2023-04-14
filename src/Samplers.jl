module Samplers
using Reexport

# Write your package code here.
include("IsingModel.jl")
include("Measurements.jl")

@reexport using .Ising
@reexport using .Measure






end
