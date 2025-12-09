using Graphs, Graphs.SimpleGraphs
using SparseArrays, LinearAlgebra
using SimpleTraits, Statistics
using Base.Threads: @threads, nthreads, threadid, Atomic, atomic_add!
using Random; Random.seed!(42)
using Dates
using CSV, DataFrames

Threads.nthreads() # Check before running

function ber_directed_divisor_graph(n::Int64, p::Float64)
    g = DiGraph(n)
    for i in 1:n, j in 2*i:i:n
        if rand() < p
            add_edge!(g, i, j)
        else
            add_edge!(g, j, i)
        end
    end
    return g
end

# PARAMETERS
n           = 2^10
ps          = 0.005:0.0001:0.015
n_replicate = 50
results_file = "00000_test_7__jl_small_2to$(Int(log2(n)))_size_lcc_and_num_scc.csv"

if !isfile(results_file)
    CSV.write(results_file, DataFrame(n=Int[], p=Float64[], rep=Int[], size_lcc=Int[], num_scc = Int[], size_plus_num_rate = Float64[]))
end

@time begin
    for p in ps
        batch_results = DataFrame(n=Int[], p=Float64[], rep=Int[], size_lcc=Int[], num_scc=Int[], size_plus_num_rate=Float64[])
        
        # Create a lock
        results_lock = ReentrantLock()
        
        Threads.@threads for rep in 1:n_replicate
            g = ber_directed_divisor_graph(n, p)

            scc_list = strongly_connected_components_tarjan(g)
            
            if !isempty(scc_list)
                # Don't build the subgraph, as the number of vertices (nv) is just the length of the list of vertices.
                scc_lengths = length.(scc_list)
                size_lcc = maximum(scc_lengths)
                
                num_scc = length(scc_list)
                size_plus_num_rate = (size_lcc + num_scc - 1) / n
                
                # Lock
                lock(results_lock) do
                    push!(batch_results, (n, p, rep, size_lcc, num_scc, size_plus_num_rate))
                end
            end
        end
        
        CSV.write(results_file, batch_results; append=true)

        @info "Completed (n,p) = ($n, $p) at $(now())"
    end
end
