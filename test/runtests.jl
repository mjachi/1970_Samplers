using Samplers
using Test
using BenchmarkTools

io = IOContext(stdout)

# Analytical solution of the critical temperature
global const IsingTc_1d = 1/(1/2*log(1+sqrt(2)))

@testset "Sanity checks" begin

    first_im = Ising.Periodic_1d(10, 1.0, 0.0)
    println(hamiltonian(first_im))
    MetropolisHastings!(first_im)
    println(hamiltonian(first_im))

    b1 = @benchmark hamiltonian(Ising.Periodic_1d(10^4, 1.0, 0.0))
    show(io, MIME("text/plain"), b1)

    b2 = @benchmark hamiltonian(Ising.Periodic_1d(10^6, 1.0, 0.0))
    show(io, MIME("text/plain"), b2)

    # second_im = IsingModel.Ising_1d(10, free, 2.0, 0.0)
    # println(hamiltonian(second_im))
end
