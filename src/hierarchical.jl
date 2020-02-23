"""
Generate hierarchical networks 
cf. Moretti & Munoz, Nat. Communi. 2013
"""

"""
    hierarchical(p::Vector{Float64}, n::Int64, seed::Int; verbose=false, write_to="")::SimpleGraph{Int}

Generate a hierarchical network with base modules of size `n` and hierarchical connection probabilities p=[p1,p2,p3...] 

This is per definition an undirected graph
# Arguments
"""
#- `log_weight(args)`: logarithmic ensemble weight function, e.g., canomical ensemble ``log\\_weight(E) = -\\beta E``
function hierarchical(M::Int64, n::Int64, seed::Int; verbose=false, write_to="")::SimpleGraph{Int}
  #...
end
