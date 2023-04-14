module Ising
"""
Metropolis-Hastings algorithm for the Ising model in 1, 2, and 3 dimensions.
"""

export Ising_1d
export periodic, free
export hamiltonian
export MetropolisHastings!

abstract type IsingModel end
abstract type Ising_1d <: IsingModel end
abstract type Ising_2d <: IsingModel end
abstract type Ising_3d <: IsingModel end

"""
  `IsingFree_1d <: Ising1d`

Struct for 1d Ising

# Fields:
-`num_spins::Int`: the number of spins in one direction

-`state::Array{Int8, 2}`: Stores the configuration over a square lattice subset.
  The array is of size `num_spins x num_spins`.

-`beta::Float64`
"""
mutable struct Free_1d <: Ising_1d
  num_spins::Int
  state::Array{Int8, 1}
  beta::T where T <: AbstractFloat
  h::S where S <: AbstractFloat
  per_site_field::Bool

  """
  Constructor for `IsingFree_1d`.

  # Arguments: 
  - `ns::Int`: Number of spins

  - `bc::BOUNDARY_CONDITIONS`: Either periodic or free

  - `beta::T`: Inverse temperature parameter

  - `h::T=0`: External field parameter
  """
  function Free_1d(
    ns::Int,
    beta::T,
    h::T=0
  ) where T <: Real
    @assert ns > 0
    new(ns, rand(Int8[-1, 1], ns), beta, h, false)
  end
end

"""
  `IsingFree_1d <: Ising1d`

Struct for 1d Ising

# Fields:
-`num_spins::Int`: the number of spins in one direction

-`state::Array{Int8, 2}`: Stores the configuration over a square lattice subset.
  The array is of size `num_spins x num_spins`.

-`beta::Float64`
"""
mutable struct Periodic_1d <: Ising_1d
  num_spins::Int
  state::Array{Int8, 1}
  beta::T where T <: AbstractFloat
  h::S where S <: AbstractFloat
  per_site_field::Bool

  """
  Constructor for `IsingPeriodic_1d`.

  # Arguments: 
  - `ns::Int`: Number of spins

  - `beta::T`: Inverse temperature parameter

  - `h::T=0`: External field parameter
  """
  function Periodic_1d(
    ns::Int,
    beta::T,
    h::T=0
  ) where T <: Real
    @assert ns > 0
    new(ns, rand(Int8[-1, 1], ns), beta, h, false)
  end
end

###########################################################################################
# Beginning of Hamiltonians

"""
  `hamiltonian`

Calculates the energy of a configuration given by `im`.

# Arguments

- `im::Ising_1d`: An `Ising_1d` instance with periodic boundary conditions

# Returns

- `Float64`: The energy of the configuration
"""
function hamiltonian(im::Periodic_1d)::Float64
  pair_sum = 0.0
  site_sum = im.h * sum(im.state)

  site_sum = sum(im.state)
  for i in 1:im.num_spins
  #@simd for i in 1:im.num_spins
    # Pair sum
    if i == im.num_spins
      pair_sum += im.state[1] * im.state[im.num_spins]
    else
      pair_sum += im.state[i] * im.state[i+1]
    end
  end

  return (-im.beta * pair_sum) - site_sum
end


"""
  `hamiltonian`

Calculates the energy of a configuration given by `im`.

# Arguments

- `im::Ising_1d`: An `Ising_1d` instance with free boundary conditions
"""
function hamiltonian(im::Free_1d)::Float64
  pair_sum = 0.0
  site_sum = im.h * sum(im.state)

  site_sum = sum(im.state)
  @simd for i in 1:im.num_spins-1
    # Pair sum
    pair_sum += im.state[i] * im.state[i+1]
  end

  return (-im.beta * pair_sum) - site_sum
end

###########################################################################################
# Beginning of sampling algs

function MetropolisHastings!(im::T, niters::Int=10^6) where {T <: Ising_1d}
  @assert niters > 0 
  niters > 100 || println("{ising.jl :: update!} -- Recommended that niters >> 1000")

  s = im.state
  ns = im.num_spins
  prob = 1/1(1+exp(-2*im.beta)*im.h)

  # TODO -- Is there a way to O(1) do this?
  oldH = hamiltonian(im)

  for i in 1:niters
    idx = rand(1:ns)
    im.state[idx] *= -1

    newH = hamiltonian(im)

    if rand() < prob || newH < oldH
      oldH = newH
    else 
      s[idx] *= -1
    end
  end
end

###########################################################################################
# Beginning of physical calculations

"""
  `boltzmann(im::Ising_1d)::Float64`

Calculates the Boltzmann measure of a 1-d configuration

# Arguments

- `im::Ising_1d`: The configuration of interest
"""
function boltzmann(im::Ising_1d)::Float64
  # TODO -- old gidas pset?

end





end