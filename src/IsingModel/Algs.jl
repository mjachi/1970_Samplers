

# Helpers
"""
Returns an array of neighbor indices
"""
function neighbors(im::T, i::Int) where T <: Ising_1d
    ns = im.num_spins
    return [CartesianIndex(ifelse(i==1, ns, i-1)), CartesianIndex(ifelse(i==ns, 1, i+1))]
end

function neighbors(im::T, i::Int) where T <: Ising_2d


end


###########################################################################################
# Metropolis
###########################################################################################






# Implementations

"""
    `Metropolis!`

Metropolis MCMC implementation for the Ising model.

# Arguments
-`im::T`: Ising_1d struct.
"""
function Metropolis!(im::T, niters::Int=10^6) where T <: Ising_1d
    if 1 > 0
        println("asdfasdfqwerqwer")
    end


end

function Metropolis!(im::T, niters::Int=10^6) where T <: Ising_2d

end

function Metropolis!(im::T, niters::Int=10^6) where T <: Ising_3d

end


###########################################################################################
# Wolff
###########################################################################################

function Wolff!(im::T, niters::Int=10^6) where T <: Ising_1d

end
function Wolff!(im::T, niters::Int=10^6) where T <: Ising_2d

end
function Wolff!(im::T, niters::Int=10^6) where T <: Ising_3d

end

###########################################################################################
# Perfect Sampling
###########################################################################################

"""
    `Perfect!(im::T) where T <: Ising_1d`

Implements the perfect sampling algorithm for the 1d Ising model.

"""
function Perfect!(im::T) where T <: Ising_1d

end
function Perfect!(im::T) where T <: Ising_2d

end
function Perfect!(im::T) where T <: Ising_3d

end

