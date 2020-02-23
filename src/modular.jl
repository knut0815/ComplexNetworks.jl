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

this is a stochastic block matrix ...
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

  #does not work it seems
  if is_directed
    #problem: need to get index k from index i
    # where k=1 for i=1:n[1], k=2 for i=n[1]+1:n[1]+n[2]
    cumsum_n = cumsum(n)
    k = 1
    for i = 1:N
      if i > cumsum_n[k]
        k += 1
      end
      l = 1
      for j = 1:N
        if j > cumsum_n[l]
          l += 1
        end
        if i!=j && rand(rng) < p[k,l]
          add_edge!(network,i,j)
        end
      end
    end
  else
    #todo
  end
  return network

  #for index in CartesianIndices(p) 
  #  i = index[1]
  #  j = index[2]
  #  if i==j
  #    max_ne = is_directed ? n[i]*(n[i]-1) : div(n[i] * (n[i] - 1), 2)
  #  else
  #    max_ne = is_directed ? n[i]*n[j]     : div(n[i]*n[j], 2)
  #  end
  #  ne = LightGraphs.SimpleGraphs.randbn(max_ne, p[index], seed=seed+i)

  #  nmin_i = i==1 ? 1 : sum(n[i-1]
  #  nmin_j = j==1 ? 1 : sum(n[j-1]
  #  
  #  nmax_i = i==lenght(n) ? N : sum(n[i+1]
  #  nmax_j = j==lenght(n) ? N : sum(n[j+1]

  #  ne_old = network.ne
  #  while network.ne - ne_old < ne
  #    source = rand(rng, nmin_i:nmax_i)
  #    dest   = rand(rng, nmin_j:nmax_j)
  #    if i==j
  #      source != dest && add_edge!(network, source, dest)
  #    else
  #      add_edge!(network, source, dest)
  #    end
  #  end
  #end
  
  
end
