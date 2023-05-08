
###########################################################################################
# Metropolis
###########################################################################################

# Helpers
"""
Returns an array of neighbor indices
"""
function neighbors(im::Periodic_1d, i::Int)
    ns = im.num_spins
    return [CartesianIndex(ifelse(i==1, ns, i-1)),
        CartesianIndex(ifelse(i==ns, 1, i+1))]
end

function neighbors(im::Periodic_2d, i::Int, j::Int)
    ns = im.num_spins
    return [CartesianIndex((ifelse(i==1, ns, i-1), j)),
        CartesianIndex((ifelse(i==ns, 1, i+1), j)),
        CartesianIndex((i, ifelse(j==1, ns, j-1))),
        CartesianIndex((i, ifelse(j==ns, 1, j+1)))]
end

function neighbors(im::Periodic_3d, i::Int,j::Int,k::Int)
    ns = im.num_spins
    return [CartesianIndex((ifelse(i==1, ns, i-1), j, k)),
        CartesianIndex((ifelse(i==ns, 1, i+1), j, k)),
        CartesianIndex((i, ifelse(j==1, ns, j-1), k)),
        CartesianIndex((i, ifelse(j==ns, 1, j+1), k)),
        CartesianIndex((i, j, ifelse(k==1, ns, k-1))),
        CartesianIndex((i, j, ifelse(k==ns, 1, k+1)))]
end

"""
To determine dE -- N.B. only need to calculate the difference based on the neighbors
to the flipped spin in Metropolis algorithm
"""

function dEnergy(im::Periodic_1d, i::Int)
    @assert i > 0
    @assert im.num_spins >= i
    nbrs = neighbors(im, i)

    bond_energy = -im.state[i]*(im.state[nbrs[1]] + im.state[nbrs[2]])
    return bond_energy
end

function dEnergy(im::Periodic_2d, i::Int,j::Int)
    @assert i > 0
    @assert im.num_spins >= i
    nbrs = neighbors(im, i, j)

    return -im.state[i,j]*(im.state[nbrs[1]] + im.state[nbrs[2]] + im.state[nbrs[3]] + im.state[nbrs[4]])
end

function dEnergy(im::Periodic_3d, i::Int,j::Int,k::Int)
    @assert i > 0
    @assert im.num_spins >= i
    nbrs = neighbors(im, i,j,k)

    return -im.state[i,j,k]*(
        im.state[nbrs[1]] + im.state[nbrs[2]] + im.state[nbrs[3]] +
        im.state[nbrs[4]] + im.state[nbrs[5]] + im.state[nbrs[6]]
    )
end

# Implementation

"""
    `Metropolis!`

Metropolis-Hastings implementation for the Ising model.

# Arguments

-`im::T`: Ising_1d struct.

-`niters::Int`: Numbers of sweeps to complete.
"""
function Metropolis!(im::T, niters::Int=1000) where T <: Ising_1d
    @assert niters > 0
    ns = im.num_spins

    @inbounds @simd for _ in 1:niters
        @inbounds @simd for _ in 1:ns
            i = rand(1:ns)

            de = -2/im.beta * dEnergy(im, i)
            if (de <= 0 || rand() < exp(-de))
                im.state[i] *= -1
                
                im.energy += de * im.beta
                im.magnetization += 2.0 * im.state[i]
            end
        end
    end
end

"""
    `Metropolis!`

Metropolis-Hastings implementation for the Ising model.

# Arguments

-`im::T`: Ising_2d struct.

-`niters::Int`: Number of sweeps to complete.
"""
function Metropolis!(im::T, niters::Int=1000) where T <: Ising_2d
    @assert niters > 0

    ns = im.num_spins

    @inbounds @simd for _ in 1:niters
        @inbounds @simd for _ in 1:ns^2
            i = rand(1:ns)
            j = rand(1:ns)

            de = -2/im.beta * dEnergy(im,i,j)
            if (de <= 0 || rand() < exp(-de))
                im.state[i,j] *= -1

                im.energy += im.beta * de
                im.magnetization += 2.0 * im.state[i,j]
            end
        end
    end
end


"""
    `Metropolis!`

Metropolis-Hastings implementation for the Ising model.

# Arguments
-`im::T`: Ising_3d struct.

-`niters::Int`: Number of sweeps to complete.

"""
function Metropolis!(im::T, niters::Int=1000) where T <: Ising_3d
    @assert niters > 0

    ns = im.num_spins
    for _ in 1:niters
        for _ in 1:ns^3
            i = rand(1:ns)
            j = rand(1:ns)
            k = rand(1:ns)
            
            de = -2/im.beta * dEnergy(im,i,j)
            if (de <= 0 || rand() < exp(-de))
                im.state[i,j,k] *= -1

                im.energy += im.beta * de
                im.magnetization += 2.0 * im.state[i,j,k]
            end
        end
    end
end


###########################################################################################
# Glauber Dynamics
###########################################################################################

"""
    `Glauber!`

Glauber dynamics implementation for the Ising model.

# Arguments

-`im::T`: Ising_1d struct.

-`niters::Int`: Numbers of sweeps to complete.
"""
function Glauber!(im::T, niters::Int=10^6) where T <: Ising_1d
    @assert niters > 0
    ns = im.num_spins

    @inbounds @simd for _ in 1:niters
        @inbounds @simd for _ in 1:ns
            i = rand(1:ns)

            de = -2/im.beta * dEnergy(im, i)
            if (rand() < (exp(-de) / (1 + exp(-de))))
                im.state[i] *= -1
                
                im.energy += de * im.beta
                im.magnetization += 2.0 * im.state[i]
            end
        end
    end
end

"""
    `Glauber!`

Glauber dynamics implementation for the Ising model.

# Arguments

-`im::T`: Ising_2d struct.

-`niters::Int`: Number of sweeps to complete.
"""
function Glauber!(im::T, niters::Int=10^6) where T <: Ising_2d
    @assert niters > 0

    ns = im.num_spins

    @inbounds @simd for _ in 1:niters
        @inbounds @simd for _ in 1:ns^2
            i = rand(1:ns)
            j = rand(1:ns)

            de = -2/im.beta * dEnergy(im,i,j)
            if (rand() < (exp(-de) / (1 + exp(-de))))
                im.state[i,j] *= -1

                im.energy += im.beta * de
                im.magnetization += 2.0 * im.state[i,j]
            end
        end
    end
end


"""
    `Glauber`

Glauber dynamics implementation for the Ising model.

# Arguments

-`im::T`: Ising_3d struct.

-`niters::Int`: Number of sweeps to complete.
"""
function Glauber!(im::T, niters::Int=10^6) where T <: Ising_3d
    @assert niters > 0

    ns = im.num_spins
    for _ in 1:niters
        for _ in 1:ns^3
            i = rand(1:ns)
            j = rand(1:ns)
            k = rand(1:ns)
            
            de = -2/im.beta * dEnergy(im,i,j)
            if (rand() < (exp(-de) / (1 + exp(-de))))
                im.state[i,j,k] *= -1

                im.energy += im.beta * de
                im.magnetization += 2.0 * im.state[i,j,k]
            end
        end
    end
end


###########################################################################################
# Wolff
###########################################################################################

# Helpers

"""
Forms the cluster to be flipped
"""
function wolff_cluster(im::T, i::CartesianIndex) where {T <: Ising_1d}
    # padd = wolff_p_add(1/im.beta)
    cluster = falses(size(im.state)...)

    queue = [i]
    while !isempty(queue)
        j = pop!(queue)
        for x in neighbors(im, j[1])
            if !(cluster[x]) && im.state[x] == im.state[j] && rand() < -expm1(-2/im.beta)
                cluster[x] = true
                push!(queue, x)
            end
        end
    end

    return cluster
end

function wolff_cluster(im::T, i::CartesianIndex) where {T <: Ising_2d}
    cluster = falses(size(im.state)...)

    queue = [i]
    while !isempty(queue)
        j = pop!(queue)
        for x in neighbors(im, j[1], j[2])
            if !(cluster[x]) && im.state[x] == im.state[j] && rand() < -expm1(-2/im.beta)
                cluster[x] = true
                push!(queue, x)
            end
        end
    end

    return cluster
end

function wolff_cluster(im::T, i::CartesianIndex) where {T <: Ising_3d}
    cluster = falses(size(im.state)...)

    queue = [i]
    while !isempty(queue)
        j = pop!(queue)
        for x in neighbors(im, j[1], j[2], j[3])
            if !(cluster[x]) && im.state[x] == im.state[j] && rand() < -expm1(-2/im.beta)
                cluster[x] = true
                push!(queue, x)
            end
        end
    end

    return cluster
end



# Implementations

"""
    `Wolff!`

Implementation of the Wolff algorithm for the 1d Ising Model
"""
function Wolff!(im::T, niters::Int=10^6) where T <: Ising_1d
    @assert niters > 0

    ns = im.num_spins

    for _ in 1:niters
        i = rand(1:ns)
        curr_coord = CartesianIndex((i))

        cluster = wolff_cluster(im, curr_coord)
        cluster_size = sum(cluster)

        im.magnetization += -2 * im.state[curr_coord] * cluster_size
        
        im.state .= ifelse.(cluster, im.state .* -1, im.state)
        im.energy = Hamiltonian(im)
    end
end

"""
    `Wolff!`

Implementation of the Wolff algorithm for the 2d Ising Model
"""
function Wolff!(im::T, niters::Int=10^6) where T <: Ising_2d
    @assert niters > 0

    ns = im.num_spins

    for _ in 1:niters
        i = rand(1:ns)
        j = rand(1:ns)
        curr_coord = CartesianIndex((i,j))

        cluster = wolff_cluster(im, curr_coord)
        cluster_size = sum(cluster)

        im.magnetization += -2 * im.state[curr_coord] * cluster_size
        
        im.state .= ifelse.(cluster, im.state .* -1, im.state)

        im.energy = Hamiltonian(im)
    end
end

"""
    `Wolff!`

Implementation of the Wolff algorithm for the 3d Ising Model
"""
function Wolff!(im::T, niters::Int=10^6) where T <: Ising_3d
    @assert niters > 0

    ns = im.num_spins

    for _ in 1:niters
        i = rand(1:ns)
        j = rand(1:ns)
        k = rand(1:ns)
        curr_coord = CartesianIndex((i,j))

        cluster = wolff_cluster(im, curr_coord)

        cluster_size = sum(cluster)

        im.magnetization += -2 * im.state[curr_coord] * cluster_size
        im.state .= ifelse.(cluster, im.state .* -1, im.state)
        im.energy = Hamiltonian(im)
    end
end

###########################################################################################
# Perfect Sampling
###########################################################################################

"""
    `Perfect!(im::T) where T <: Ising_1d`

Implements the perfect sampling algorithm for the 1d Ising model.

"""
function Perfect!(im::T) where T <: Ising_1d
    @assert im.beta > 0

end
function Perfect!(im::T) where T <: Ising_2d
    @assert im.beta > 0

end
function Perfect!(im::T) where T <: Ising_3d
    @assert im.beta > 0

end

###########################################################################################
# Coupling from the past (CFTP)
###########################################################################################

function CFTP(im::T) where T <: Ising_1d
    @assert im.beta > 0

end

function CFTP(im::T) where T <: Ising_2d
    @assert im.beta > 0

end

function CFTP(im::T) where T <: Ising_3d
    @assert im.beta > 0

end