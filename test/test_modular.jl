
using LinearAlgebra
using Random
using Statistics
using LightGraphs

#TODO: reexport LightGraphs functions for the use of the module
#
#TODO: test with modules of size 1 and heterogeneous (random?) p matrix
#TODO: test with modules of different sizes

#takes ~3s on laptop
function test_modular_uniform(;verbose=false)
    pass = true

    m = 10  #number of modules
    p_in  = 0.5
    p_out = 0.1
    p  = p_out*(ones(m,m)-I) + I*p_in
    n0 = 50

    N = n0*m
    module_ranges = []
    for i in 1:m
      push!(module_ranges, (i-1)*n0+1 : i*n0) 
    end

    seed_range = 1000:1299

    ###########################################################################
    num_p_in_within_1sigma = 0
    num_p_out_within_1sigma = 0
    num_p_in_larger = 0
    num_p_out_larger = 0
    for seed in seed_range
      g = ComplexNetworks.modular_uniform(p,n0,is_directed=true, seed=seed) 
      num_connections_within_module  = zeros(Float64, N)
      num_connections_outside_module = zeros(Float64, N)
      for range in module_ranges
        for i = range
          for j in outneighbors(g,i)
            if j in range
              num_connections_within_module[i] += 1
            else
              num_connections_outside_module[i] += 1
            end
          end
        end
      end
      mean_num_connections_within_module = mean(num_connections_within_module)
      std_num_connections_within_module = std(num_connections_within_module, mean=mean_num_connections_within_module)
      mean_num_connections_outside_module = mean(num_connections_outside_module)
      std_num_connections_outside_module = std(num_connections_outside_module, mean=mean_num_connections_outside_module)
    
      #normalization depends on number of possible connection partners (self-connections are categorically excluded)
      p_in_meas  = mean_num_connections_within_module  / (n0-1) 
      p_in_err   = std_num_connections_within_module   / (n0-1) / sqrt(N-1)
      p_out_meas = mean_num_connections_outside_module / (n0*(m-1)) 
      p_out_err  = std_num_connections_outside_module  / (n0*(m-1)) / sqrt(N-1) 

      num_p_in_larger += p_in_meas > p_in ? 1 : 0
      num_p_out_larger += p_out_meas > p_out ? 1 : 0

      num_p_in_within_1sigma += abs(p_in_meas - p_in) < p_in_err ? 1 : 0
      num_p_out_within_1sigma += abs(p_out_meas - p_out) < p_out_err ? 1 : 0
    end

    percent_p_in_within_1sigma = (num_p_in_within_1sigma/length(seed_range)*100)
    percent_p_out_within_1sigma = (num_p_out_within_1sigma/length(seed_range)*100)
    percent_p_in_larger = (num_p_in_larger/length(seed_range)*100)
    percent_p_out_larger = (num_p_out_larger/length(seed_range)*100)
    
    #statistically 1 sigma should include 68.3%->68%(floor) and 50% should be larger
    result  = floor(percent_p_in_within_1sigma)  in 66:70 && floor(percent_p_in_larger)  in 47:53 
    result &= floor(percent_p_out_within_1sigma) in 66:70 && floor(percent_p_out_larger) in 47:53 

    if verbose 
      println("... directed, all $(m) modules have same number of nodes [$(result)]")
      println("...   percent p_in within 1sigma  = $(percent_p_in_within_1sigma)") 
      println("...   percent p_out within 1sigma = $(percent_p_out_within_1sigma)") 
      println("...   percent p_in larger         = $(percent_p_in_larger)") 
      println("...   percent p_out larger        = $(percent_p_out_larger)") 
    end

    pass &= result

    ###########################################################################
    ##TODO: abstract into function
    num_p_in_within_1sigma = 0
    num_p_out_within_1sigma = 0
    num_p_in_larger = 0
    num_p_out_larger = 0
    for seed in seed_range
      g = ComplexNetworks.modular_uniform(p,n0,is_directed=false, seed=seed) 
      num_connections_within_module  = zeros(Float64, m)
      num_connections_outside_module = zeros(Float64, m)
      for (mi,range) in enumerate(module_ranges)
        for i = range
          for j in outneighbors(g,i)
            mj = mi==m ? 1 : mi+1
            if j in range
              num_connections_within_module[mi] += 1
            elseif j in module_ranges[mj]
              num_connections_outside_module[mi] += 1
            end
          end
        end
      end
      mean_num_connections_within_module = mean(num_connections_within_module)
      std_num_connections_within_module = std(num_connections_within_module, mean=mean_num_connections_within_module)
      mean_num_connections_outside_module = mean(num_connections_outside_module)
      std_num_connections_outside_module = std(num_connections_outside_module, mean=mean_num_connections_outside_module)
    
      #normalization depends on number of possible connection partners (self-connections are categorically excluded)
      p_in_meas  = mean_num_connections_within_module  / (n0*(n0-1)) 
      p_in_err   = std_num_connections_within_module   / (n0*(n0-1)) / sqrt(m-1)
      p_out_meas = mean_num_connections_outside_module / (n0*n0) 
      p_out_err  = std_num_connections_outside_module  / (n0*n0) / sqrt(m-1) 

      num_p_in_larger += p_in_meas > p_in ? 1 : 0
      num_p_out_larger += p_out_meas > p_out ? 1 : 0

      num_p_in_within_1sigma += abs(p_in_meas - p_in) < p_in_err ? 1 : 0
      num_p_out_within_1sigma += abs(p_out_meas - p_out) < p_out_err ? 1 : 0
    end

    percent_p_in_within_1sigma = (num_p_in_within_1sigma/length(seed_range)*100)
    percent_p_out_within_1sigma = (num_p_out_within_1sigma/length(seed_range)*100)
    percent_p_in_larger = (num_p_in_larger/length(seed_range)*100)
    percent_p_out_larger = (num_p_out_larger/length(seed_range)*100)
    
    #statistically 1 sigma should include 68.3%->68%(floor) and 50% should be larger
    result  = floor(percent_p_in_within_1sigma) in 66:70 && floor(percent_p_in_larger) in 48:52 
    result &= floor(percent_p_out_within_1sigma) in 66:70 && floor(percent_p_out_larger) in 48:52 

    if verbose 
      println("... undirected, all $(m) modules have same number of nodes [$(result)]")
      println("...   percent p_in within 1sigma  = $(percent_p_in_within_1sigma)") 
      println("...   percent p_out within 1sigma = $(percent_p_out_within_1sigma)") 
      println("...   percent p_in larger         = $(percent_p_in_larger)") 
      println("...   percent p_out larger        = $(percent_p_out_larger)") 
    end
    pass &= result

    return pass
end
