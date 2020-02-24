"""
Generate modular networks 

Literature:
Lancichinetti, A., Fortunato, S. & Radicchi, F. Benchmark graphs for testing community detection algorithms. Phys. Rev. E 78,
046110, https://doi.org/10.1103/PhysRevE.78.046110 (2008).

Andrea Lancichinetti and Santo Fortunato
Benchmarks for testing community detection algorithms on directed and weighted graphs with overlapping communities
Phys. Rev. E 80, 016118 – Published 31 July 2009
"""

"""
    modular(p::Array{Float64,2}, n::Int64, seed::Int; verbose=false, write_to="")::SimpleGraph{Int}

Generate a modular network of base modules of size `n` with connection probability matrix p 

This is per definition an undirected graph
# Arguments
"""
function modular(P_degree, P_size, m::Int64, n::Int64, seed::Int; verbose=false, write_to="")::SimpleGraph{Int}
end


"""
network with modules of sizes `n` with edges drawn with probability `p[i,j]`
"""
function modular_uniform(p::Array{Float64,2}, n::Int64; seed::Integer=-1, is_directed=false)
  list_n = ones(Integer, size(p)[1])*n
  modular_uniform(p, list_n, seed=seed, is_directed=is_directed)
end

"""
network with modules of sizes `n[k]` with edges drawn with probability `p[k,l]`
excludes self-connections 
example:
using LinearAlgebra
p=0.1*(ones(3,3)-I) + I*0.5
g=modular_uniform(p, 10, is_directed=true)
am = Matrix(adjacency_matrix(g))

this is similar to stochastic block matrix ...
"""
function modular_uniform(p::Array{Float64,2}, n::Vector{Int64}; seed::Integer=-1, is_directed=false)
  @assert size(p)[1] == size(p)[2]  #square matrix
  @assert size(p)[1] == length(n)

  N = sum(n)
  network = is_directed ? SimpleDiGraph(N) : SimpleGraph(N)

  if seed >= 0
    rng = MersenneTwister(seed)
  else
    rng = Random.GLOBAL_RNG
  end

  #we need the cumsum in order to evaluate to which module each index belongs
  cumsum_n = cumsum(n)
  module_i = 1
  for i = 1:N
    module_i += i <= cumsum_n[module_i] ? 0 : 1
    #this is a sufficient start because even though j_min can be i+1, module_j
    #is checked at beginning of loop and would increase in case of violation
    module_j = module_i
    for j = i+1:N
      module_j += j <= cumsum_n[module_j] ? 0 : 1
      if rand(rng) < p[module_i,module_j]
        add_edge!(network,i,j)
      end
      if is_directed == true # in case directed also draw reverse edge with probability
        if rand(rng) < p[module_j,module_i]
          add_edge!(network,j,i)
        end
      end
    end
  end
  return network
end

