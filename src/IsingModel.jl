module Ising
using .Threads

"""
Metropolis-Hastings algorithm for the Ising model in 1, 2, and 3 dimensions.
"""

export Ising_1d
export periodic, free
export hamiltonian
export total_magnetization
export MetropolisHastings!

abstract type IsingModel end
abstract type Ising_1d <: IsingModel end
abstract type Ising_2d <: IsingModel end
abstract type Ising_3d <: IsingModel end

global const crit_temp_2d = (2 / log(1 + sqrt(2)));

##########################################################################################
# Configuration struct definitions

"""
  `Free_1d <: Ising1d`

Struct for 1d Ising with free boundary conditions

# Fields:
-`num_spins::Int`: the number of spins in one direction

-`state::Array{Int8, 1}`: Stores the configuration over a square lattice subset.
  The array is of size `num_spins`.

-`beta::T`: Inverse temperature parameter

-`h::S`: External magnetic field parameter; default 0 in constructor
"""
mutable struct Free_1d <: Ising_1d
  num_spins::Int
  state::Array{Int8, 1}
  beta::T where T <: AbstractFloat
  h::S where S <: AbstractFloat

  """
  Constructor for `Free_1d`.

  # Arguments: 
  - `ns::Int`: Number of spins

  - `beta::T`: Inverse temperature parameter

  - `h::T=0`: External magnetic field parameter; default 0
  """
  function Free_1d(
    ns::Int,
    beta::T,
    h::T=0
  ) where T <: Real
    @assert ns > 0
    new(ns, rand(Int8[-1, 1], ns), beta, h)
  end
end

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
  `Free_2d <: Ising2d`

Struct for 2d Ising

# Fields:
-`num_spins::Int`: the number of spins in one direction

-`state::Array{Int8, 2}`: Stores the configuration over a square lattice subset.
  The array is of size `num_spins x num_spins`.

-`beta::T`: Inverse temperature parameter

-`h::S`: External magnetic field parameter
"""
mutable struct Free_2d <: Ising_2d
  num_spins::Int
  state::Array{Int8, 2}
  beta::T where T <: AbstractFloat
  h::S where S <: AbstractFloat

  """
  Constructor for `Free_2d`.

  # Arguments: 

  - `ns::Int`: Number of spins

  - `beta::T`: Inverse temperature parameter

  - `h::T=0`: External field parameter
  """
  function Free_2d(
    ns::Int,
    beta::T,
    h::T=0
  ) where T <: Real
    @assert ns > 0
    new(ns, rand(Int8[-1, 1], (ns, ns)), beta, h)
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
  `Free_3d <: Ising3d`

Struct for 3d Ising

# Fields:
-`num_spins::Int`: the number of spins in one direction

-`state::Array{Int8, 3}`: Stores the configuration over a square lattice subset.
  The array is of size `num_spins x num_spins x num_spins`.

-`beta::T`: Inverse temperature parameter

-`h::S`: External magnetic field parameter
"""
mutable struct Free_3d <: Ising_3d
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
  function Free_3d(
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

  @inbounds @simd for i in 1:im.num_spins
    pair_sum += im.state[i] * im.state[i == im.num_spins ? 1 : i+1]
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

  @inbounds @simd for i in 1:im.num_spins-1
    # Pair sum
    pair_sum += im.state[i] * im.state[i+1]
  end

  return (-im.beta * pair_sum) - site_sum
end

"""
Free 2d hamiltonian
"""
function hamiltonian(im::Free_2d)::Float64
  pair_sum = 0.0
  site_sum = im.h * sum(im.state)

  for i in 1:im.num_spins
    for j in 1:im.num_spins
      pair_sum += i+1 > im.num_spins ? 0 : im.state[i,j] * im.state[i+1, j]
      pair_sum += j+1 > im.num_spins ? 0 : im.state[i,j] * im.state[i, j+1]
    end
  end

  return (-im.beta * pair_sum) - site_sum
end

"""
Periodic 2d hamiltonian
"""
function hamiltonian(im::Periodic_2d)::Float64
  pair_sum = 0.0
  site_sum = im.h * sum(im.state)

  for i in 1:im.num_spins
    for j in 1:im.num_spins
      pair_sum += im.state[i,j] * (i == im.num_spins ? im.state[1,j] : im.state[i+1, j])
      pair_sum += im.state[i,j] * (j == im.num_spins ? im.state[i,1] : im.state[i, j+1])
    end
  end

  return (-im.beta * pair_sum) - site_sum
end

##########################################################################################

function total_magnetization(im::T) where T <: IsingModel
  return mean(im.state)
end

function neighbors(im::Periodic_1d, i::Int)
  ns = im.num_spins
  return [CartesianIndex(ifelse(i==1, ns, i-1)), CartesianIndex(ifelse(i==ns, 1, i+1))]
end

function neighbors(im::Periodic_2d, i,j)
  ns = im.num_spins
  return [CartesianIndex((ifelse(i==1, ns, i-1), j)),
    CartesianIndex((ifelse(i==ns, 1, i+1), j)),
    CartesianIndex((i, ifelse(j==1, ns, j-1))),
    CartesianIndex((i, ifelse(j==ns, 1, j+1)))]
end

function neighbors(im::Periodic_3d, i,j,k)
  return [CartesianIndex((ifelse(i==1, ns, i-1), j, k)),
    CartesianIndex((ifelse(i==ns, 1, i+1), j, k)),
    CartesianIndex((i, ifelse(j==1, ns, j-1), k)),
    CartesianIndex((i, ifelse(j==ns, 1, j+1), k)),
    CartesianIndex((i, j, ifelse(k==1, ns, k-1))),
    CartesianIndex((i, j, ifelse(k==ns, 1, k+1)))]
end

###########################################################################################
# Beginning of sampling algs

"""
Metropolis-Hastings algorithm for 1d

ASSUMES PERIODIC, 0 FIELD
"""
function MetropolisHastings!(im::T, niters::Int=10^6) where {T <: Ising_1d}
  @assert niters > 0 

  ns = im.num_spins

  for _ in 1:niters
    for _ in 1:ns
      i = rand(1:ns)

      energy = bond_energ
    end
  end

end

function bond_energy(im::T, i, j) where T <: Ising_2d
  nbrs = neighbors(im, i,j)
  
  return -im.state[i,j]*(im.state[nbrs[1]] + im.state[nbrs[2]] + im.state[nbrs[3]] + im.state[nbrs[4]])
end


"""
Metropolis-Hastings algorithm for 2d

ASSUMES PERIODIC, ASSUME 0 FIELD
"""
function MetropolisHastings!(im::T, niters::Int=10^6) where {T <: Ising_2d}
  @assert niters > 0

  ns = im.num_spins

  for _ in 1:niters
    for _ in 1:ns^2
      i = rand(1:ns)
      j = rand(1:ns)
  
      energy = bond_energy(im, i,j) / im.beta
      de = -2 * (energy + 0)

      if (de <= 0 || rand() < exp(-de / crit_temp))
        im.state[i,j] *= -1
        im.energy += de / (ns^2)
        im.magnetization += 2.0 * im.state[i,j] * im.state[ns, ns]
      end
    end
  end

end

"""
Metropolis-Hastings algorithm for 3d
"""
function MetropolisHastings!(im::T, niters::Int=10^6) where {T <: Ising_2d}
  @assert niters > 0

end

"""
Wolff clustering for 1d
"""
function Wolff!(im::T, niters::Int=10^6) where {T <: Ising_1d}

end

"""
Wolff clustering
"""
function Wolff!(im::T, niters::Int=10^6) where {T <: Ising_2d}

end

"""

"""
function Wolff!(im::T, niters::Int=10^6) where {T <: Ising_3d}

end

###########################################################################################
# Beginning of physical calculations




end
