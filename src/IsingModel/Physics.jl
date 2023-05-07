"""
Implements several physical calculations
"""

module Physics

export SpecificHeat, Susceptibility

"""
    `SpecificHeat`

Equilibrates and calculates the magnetic susceptibility for an Ising model struct
"""
function SpecificHeat(
    im::T;
    algorithm::Symbol=:Metropolis,
    nsweeps::Int=10^6,
    ntherm::Int=10^5,
    measurement_interval::Int=10^2
) where T <: IsingModel
    E_acc = Measurements.Accumulator()   
    E2_acc = Measurements.Accumulator()

    # Equilibration steps
    


end


"""
    `Susceptibility(im::T) where T <: IsingModel`

Equilibrates and calculates the magnetic susceptibility for an Ising model struct.
"""
function Susceptibility(
    im::T,
    nsweeps::Int=10^6,
    nterhm::Int=10^5,
    measurement_interval::Int=10^2
) where T <: IsingModel


end

end