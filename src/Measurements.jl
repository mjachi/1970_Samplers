module Measurements

export Accumulator, add!, mean
"""
  `Accumulator`

Struct for accumulating an observable and count for taking the empirical mean.
"""
mutable struct Accumulator
  count::Dict{String, UInt64}
  data::Dict{String, Any}

  function Accumulator()
    new(Dict{String, UInt64}(), Dict{String, Any}())
  end
end

"""
  `add!`

Add an observation.
"""
function add!(acc::Accumulator, name::String, data::T) where T <: AbstractFloat
  if haskey(acc.count, name)
    acc.count[name] += 1
    acc.data[name] += data
  else
    acc.count[name] = 1
    acc.data[name] = copy(data)
  end
end

"""
  `mean`

Take the empirical mean according to the struct
"""
function mean(acc::Accumulator, name::String)::Float64
  return acc.data[name] / acc.count[name]
end

end