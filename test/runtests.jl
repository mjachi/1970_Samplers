using Samplers
using Test
using BenchmarkTools
using Plots

io = IOContext(stdout)

# Analytical solution of the critical temperature
global const IsingTc_1d = 2/(log(1+sqrt(2)))

@testset "Sanity checks" begin

  


end
