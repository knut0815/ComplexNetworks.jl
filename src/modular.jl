"""
Generate modular networks 

Literature:
Lancichinetti, A., Fortunato, S. & Radicchi, F. Benchmark graphs for testing community detection algorithms. Phys. Rev. E 78,
046110, https://doi.org/10.1103/PhysRevE.78.046110 (2008).

Andrea Lancichinetti and Santo Fortunato
Benchmarks for testing community detection algorithms on directed and weighted graphs with overlapping communities
Phys. Rev. E 80, 016118 – Published 31 July 2009
"""
using LightGraphs
using Random

"""
    modular(p::Array{Float64,2}, n::Int64, seed::Int; verbose=false, write_to="")::SimpleGraph{Int}

Generate a modular network of base modules of size `n` with connection probability matrix p 

This is per definition an undirected graph
# Arguments
"""
function modular(m::Int64, n::Int64, seed::Int; verbose=false, write_to="")::SimpleGraph{Int}
end
