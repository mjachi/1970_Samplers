using Samplers
using Test
using BenchmarkTools
using Plots

io = IOContext(stdout)

@testset "1d tests" begin
    include("OneDim.jl")
end

@testset "2d tests" begin
    include("TwoDim.jl")
end