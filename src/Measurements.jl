module Measurements
"""
  `Accumulator`

Struct for accumulating an observable and count for taking the empirical mean.
"""
mutable struct Accumulator{T}
  count::UInt64
  data::T

  function Accumulator(initial::T) where {T <: AbstractFloat}
    new{T}(UInt64(1), initial)
  end
end

"""
  `add!`

Add an observation.
"""
function add!(acc::Accumulator{T}, data::T) where T <: AbstractFloat
  acc.count += 1
  acc.data .+= data
end

"""
  `mean`

Take the empirical mean according to the struct
"""
function mean(acc::Accumulator{T}) where T <: AbstractFloat
  return acc.data / acc.count
end

end