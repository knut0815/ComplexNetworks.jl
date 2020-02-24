using ComplexNetworks
using Test

@testset "ComplexNetworks.jl" begin
    # Write your own tests here.
    include("test_modular.jl")
    @test test_stochastic_block_model()
end
