using Samplers
using Test
using BenchmarkTools

@testset "Sanity checks" begin
    # Write your tests here.

    asdf = IsingModel.Ising_1d(10, periodic, 2.0, 0.0)
end
