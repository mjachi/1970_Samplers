module Samplers
using Reexport

# Write your package code here.
include("ising.jl")

@reexport using .IsingModel

"""
Defines measurement structures to Julian-ly keep track of observables
"""


end
