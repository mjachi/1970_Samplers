module IsingModel
"""
Metropolis-Hastings algorithm for the Ising model in 1, 2, and 3 dimensions.
"""

export Ising_1d
export periodic, free



abstract type IsingMetropolis end
@enum BOUNDARY_CONDITIONS begin
  periodic=1 
  free=2
end



"""
  `IsingMetropolis_1d <: IsingMetropolis`

Struct for 1d Ising

# Fields:
-`num_spins::Int`: the number of spins in one direction

-`state::Array{Int8, 2}`: Stores the configuration over a square lattice subset.
  The array is of size `num_spins x num_spins`.

-`beta::Float64`
"""
mutable struct Ising_1d <: IsingMetropolis
  num_spins::Int
  state::Array{Int8, 1}
  beta::T where T <: AbstractFloat
  J::Array{Int8, 1}
  bc::BOUNDARY_CONDITIONS

  per_site_field::Bool

  """
  asdfasdf
  """
  function Ising_1d(
    ns::Int,
    bc::BOUNDARY_CONDITIONS,
    beta::T, J::Array{Int8, 1} = zeros(ns)
  ) where T <: Real
    @assert ns > 0
    new(ns, rand(Int8[-1, 1], ns), beta, J, bc, true)
  end

  """
  asdfasdf
  """
  function Ising_1d(ns::Int, bc::BOUNDARY_CONDITIONS, beta::T, h::T=0) where T <: Real
    @assert ns > 0
    new(ns, rand(Int8[-1, 1], ns), beta, h * ones(ns), bc, false)
  end
end

###########################################################################################

"""
  `hamiltonian`

"""
function hamiltonian(im::Ising_1d)
  if im.bc == periodic
    return hamiltonian_periodic(im)
  end
end

"""
  `hamiltonian_periodic`

Calculates the energy of a configuration given by `im`, assuming
`im` has periodic boundary conditions.

# Arguments

- `im::Ising_1d`: An `Ising_1d` instance with periodic boundary conditions
"""
function hamiltonian_periodic(im::Ising_1d)::Float64
  pair_sum = 0.0
  site_sum = 0.0

  if !im.per_site_field
    site_sum = im.J[1] * sum(im.state)
  end

  site_sum = sum(im.state)
  for i in 1:im.num_spins
    # Pair sum
    if i == im.num_spins
      pair_sum += im.state[1] * im.state[im.num_spins]
    else
      pair_sum += im.state[i] * im.state[i+1]
    end

    # site sum
    if im.per_site_field
      site_sum += im.J[i] * im.state[i]
    end
  end

  return (-s.beta * pair_sum) - site_sum
end


"""
  `hamiltonian_free`

Calculates the energy of a configuration given by `im`, assuming
`im` has free boundary conditions.

# Arguments

- `im::Ising_1d`: An `Ising_1d` instance with free boundary conditions
"""
function hamiltonian_free(im::Ising_1d)::Float64
  pair_sum = 0.0
  site_sum = 0.0

  if !im.per_site_field
    site_sum = im.J[1] * sum(im.state)
  end

  site_sum = sum(im.state)
  for i in 1:im.num_spins-1
    # Pair sum
    pair_sum += im.state[i] * im.state[i+1]

    # site sum
    if im.per_site_field
      site_sum += im.J[i] * im.state[i]
    end
  end

  return (-s.beta * pair_sum) - site_sum
end



###########################################################################################





end


