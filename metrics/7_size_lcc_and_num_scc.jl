# size_lcc.jl

using Graphs, Graphs.SimpleGraphs
using SparseArrays, LinearAlgebra
using SimpleTraits, Statistics
using Base.Threads: @threads, nthreads, threadid, Atomic, atomic_add!
using Random; Random.seed!(42)
using Dates
using CSV, DataFrames

Threads.nthreads() # Check number of threads!

# function for graph generation
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
n           = 2^20
ps          = 0.29:0.01:0.5
n_replicate = 1_00

results_file = "7_2to$(Int(log2(n)))_size_lcc_and_num_scc.csv"

# Write header only once at the beginning
if !isfile(results_file)
    CSV.write(results_file, DataFrame(n=Int[], p=Float64[], rep=Int[], size_lcc=Int[], num_scc = Int[], size_plus_num_rate = Int[]))
end

@time begin
    for p in ps
        # Create separate DataFrame for each thread
        results_per_thread = [DataFrame(n=Int[], p=Float64[], rep=Int[], size_lcc=Int[], num_scc = Int[], size_plus_num_rate = Float64[]) for _ in 1:Threads.nthreads()]
        
        # Parallelize the loop
        Threads.@threads for rep in 1:n_replicate
            # ID of current thread
            tid = Threads.threadid()
            
            g = ber_directed_divisor_graph(n, p)

            # find largest SCC
            scc_list   = strongly_connected_components_tarjan(g)
            comp_sizes = length.(scc_list)
            verts      = scc_list[argmax(comp_sizes)]
            subg, _    = induced_subgraph(g, verts)
            size_lcc   = nv(subg) # number of vertices of the largest SCC
            num_scc    = length(scc_list)
            size_plus_num_rate = (size_lcc + num_scc - 1) / n # the "-1" comes in!! If LCC=2^n, num scc will be 1 so total 2^10+1. 
            
            # Push the result to the DataFrame specific to the thread
            push!(results_per_thread[tid], (n, p, rep, size_lcc, num_scc, size_plus_num_rate))
        end
        
        # 3. After the loop, merge the results from all threads.
        temp_results = vcat(results_per_thread...)
        
        # Write the combined results for this p to the file.
        CSV.write(results_file, temp_results; append=true)

        @info "Completed (n,p) = ($n, $p) at $(now())"
    end
end








