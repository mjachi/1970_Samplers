using Samplers
###########################################################################################
# Analytical values
###########################################################################################

global const CriticalTemperature_2d = 2 / (log(1 + sqrt(2)))


"""

"""
function onsager_magnetization(beta::T) where {T <: AbstractFloat}
    @assert beta > 0
    return max(1-csch(2*beta)^4, 0)^(1/oftype(beta, 8))
end


"""

"""
function onsager_internal_energy(beta::T) where {T <: AbstractFloat}
    k = 2 * tanh(2 * beta) / cosh(2 * beta)
    j = 2 * tanh(2 * beta) - 1
    K = SpecialFunctions.ellipk(k^2)

    return -coth(2 * beta) * (1 + 2/π * j * K)
end

"""

"""
function onsager_heat_capacity(beta::T) where {T <: AbstractFloat}
    k = 2 * tanh(2 * beta) / cosh(2 * beta)
    K = SpecialFunctions.ellipk(k^2)
    E = SpecialFunctions.ellipe(k^2)
    j = 2 * tanh(2 * beta)^2 - 1
    return beta^2 * coth(2*beta)^2 * (2/π) * (((j - 1//2)^2 + 7//4) * K - 2E - (1 - j) * π / 2)
end

##########################################################################################
# Tests
##########################################################################################
@testset "Average magnetization" begin
    im1 = Periodic_2d(250, 100.0)
    im2 = Periodic_2d(250, 10.0)
    im3 = Periodic_3d(250, 5.0)
    Metropolis!(im1)
    Metropolis!(im2)
    Metropolis!(im3)
    @test abs(im1.average_magnetization) <= 0.001
    @test abs(im2.average_magnetization) <= 0.001
    @test abs(im3.average_magnetization) <= 0.001
end

@testset "Susceptibility" begin

end

