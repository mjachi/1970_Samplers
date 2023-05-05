module Ising

export IsingMode, Ising_1d, Ising_2d, Ising_3d
export Metropolis!, Wolff!, Perfect!

abstract type IsingModel end
abstract type Ising_1d <: IsingModel end
abstract type Ising_2d <: IsingModel end
abstract type Ising_3d <: IsingModel end


"""
  `Periodic_1d <: Ising1d`

Struct for 1d Ising with periodic boundary conditions

# Fields:
-`num_spins::Int`: the number of spins in one direction

-`state::Array{Int8, 1}`: Stores the configuration over a square lattice subset.
  The array is of size `num_spins`.

-`beta::T`: Inverse temperature parameter

-`h::S`: External magnetic field parameter; default 0 in constructor
"""
mutable struct Periodic_1d <: Ising_1d
  num_spins::Int
  state::Array{Int8, 1}
  beta::T where T <: AbstractFloat
  h::S where S <: AbstractFloat

  """
  Constructor for `IsingPeriodic_1d`.

  # Arguments: 
  - `ns::Int`: Number of spins

  - `beta::T`: Inverse temperature parameter

  - `h::T=0`: External magnetic field parameter
  """
  function Periodic_1d(
    ns::Int,
    beta::T,
    h::T=0
  ) where T <: Real
    @assert ns > 0
    new(ns, rand(Int8[-1, 1], ns), beta, h)
  end
end

"""
  `Periodic_2d <: Ising2d`

Struct for 2d Ising

# Fields:
-`num_spins::Int`: the number of spins in one direction

-`state::Array{Int8, 2}`: Stores the configuration over a square lattice subset.
  The array is of size `num_spins x num_spins`.

-`beta::T`: Inverse temperature parameter

-`h::S`: External magnetic field parameter
"""
mutable struct Periodic_2d <: Ising_2d
  num_spins::Int
  state::Array{Int8, 2}
  beta::T where T <: AbstractFloat
  h::S where S <: AbstractFloat

  """
  Constructor for `Periodic_2d`.

  # Arguments: 
  - `ns::Int`: Number of spins

  - `beta::T`: Inverse temperature parameter

  - `h::T=0`: External field parameter
  """
  function Periodic_2d(
    ns::Int,
    beta::T,
    h::T=0
  ) where T <: Real
    @assert ns > 0
    new(ns, rand(Int8[-1, 1], (ns, ns)), beta, h)
  end
end

"""
  `Periodic_3d <: Ising3d`

Struct for 3d Ising

# Fields:
-`num_spins::Int`: the number of spins in one direction

-`state::Array{Int8, 3}`: Stores the configuration over a square lattice subset.
  The array is of size `num_spins x num_spins x num_spins`.

-`beta::T`: Inverse temperature parameter

-`h::S`: External magnetic field parameter
"""
mutable struct Periodic_3d <: Ising_3d
  num_spins::Int
  state::Array{Int8, 3}
  beta::T where T <: AbstractFloat
  h::S where S <: AbstractFloat

  """
  Constructor for `Periodic_3d`.

  # Arguments: 
  - `ns::Int`: Number of spins

  - `beta::T`: Inverse temperature parameter

  - `h::T=0`: External field parameter
  """
  function Periodic_3d(
    ns::Int,
    beta::T,
    h::T=0
  ) where T <: Real
    @assert ns > 0
    new(ns, rand(Int8[-1, 1], (ns, ns)), beta, h)
  end
end

include("Physics.jl")
include("Algs.jl")
include("Analysis.jl")

end