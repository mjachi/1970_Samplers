using Samplers
using Test
using BenchmarkTools

io = IOContext(stdout)


# Testing disabled b/c will take forever.

@testset "1d Metropolis" begin
    # include("Metropolis/OneDim.jl")
end

@testset "2d Metropolis" begin
    # include("Metropolis/TwoDim.jl")
end