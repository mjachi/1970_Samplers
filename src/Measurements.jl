module Measure
"""
Implements QoL structs with which to keep track of observables.
"""


"""
  `Measurements`

Structure for retaining information over time.
"""
struct Measurements
    configurations::Array{Int8, 3}
    energies::Array{Float64, 1}
    magnetization::Array{Float64, 1}
end





end