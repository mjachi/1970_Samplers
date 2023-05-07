module Ising

export IsingMode, Ising_1d, Ising_2d, Ising_3d
export Periodic_1d, Periodic_2d, Periodic_3d
export Metropolis!, Wolff!, Perfect!

abstract type IsingModel end
abstract type Ising_1d <: IsingModel end
abstract type Ising_2d <: IsingModel end
abstract type Ising_3d <: IsingModel end

"""
  `Periodic_1d <: Ising1d`

Struct for 1d Ising with periodic boundary conditions
"""
mutable struct Periodic_1d <: Ising_1d
  num_spins::Int
  state::Array{Int8, 1}
  beta::T where T <: AbstractFloat

  energy::S where S <: AbstractFloat
  magnetization::U where U <: AbstractFloat

  """
  Constructor for `IsingPeriodic_1d`.
  """
  function Periodic_1d(
    ns::Int,
    beta::T,
  ) where T <: AbstractFloat
    @assert ns > 0
    @assert beta > 0

    init = rand(Int8[-1,1], ns)

    energy = 0.0
    @inbounds @simd for i in 1:ns
      energy += init[i] * init[i == ns ? 1 : i+1]
    end

    new(ns, init, beta, beta*energy, sum(init))
  end
end

"""
  `Periodic_2d <: Ising_2d`

Struct for 2d Ising model with periodic boundary conditions
"""
mutable struct Periodic_2d <: Ising_2d
  num_spins::Int
  state::Array{Int8, 2}
  beta::T where T <: AbstractFloat

  energy::S where S <: AbstractFloat
  magnetization::U where U <: AbstractFloat

  """
  Constructor for `Periodic_2d`.
  """
  function Periodic_2d(
    ns::Int,
    beta::T,
  ) where T <: Real
    @assert ns > 0
    @assert beta > 0

    init = rand(Int8[-1,1], ns)

    energy = 0.0
    # N.B. Due to periodicity, it is sufficient to only look "north-east"
    @inbounds @simd for i in 1:ns
      @inbounds @simd for j in 1:ns
        energy += init[i,j] * (init[i==ns ? 1 : i+1, j] + init[i, j==ns ? 1 : j+1])
      end
    end

    mag = sum(init)

    new(ns, init, beta, beta*energy, mag)
  end
end

"""
  `Periodic_3d <: Ising3d`

Struct for 3d Ising model with periodic boundary conditions
"""
mutable struct Periodic_3d <: Ising_3d
  num_spins::Int
  state::Array{Int8, 3}
  beta::T where T <: AbstractFloat

  energy::S where S <: AbstractFloat
  magnetization::U where U <: AbstractFloat

  """
  Constructor for `Periodic_3d`.
  """
  function Periodic_3d(
    ns::Int,
    beta::T,
  ) where T <: Real
    @assert ns > 0
    @assert beta > 0

    init = rand(Int8[-1,1], (ns, ns, ns))

    energy = 0.0
    # N.B. Due to periodicity, it is sufficient to only look "north-east"
    @inbounds @simd for i in 1:ns
      @inbounds @simd for i in 1:ns
        @inbounds @simd for i in 1:ns
          nbr_sum = init[i==ns ? 1 : i+1,j,k] +
            init[i,j==ns ? 1 : j+1,k] + init[i,j,k==ns ? 1 : k+1]

          energy += init[i,j] * nbr_sum
        end
      end
    end

    mag = sum(init)

    new(ns, init, beta, beta*energy, mag)
  end
end

include("Algs.jl")
include("Analysis.jl")
include("Physics.jl")

end