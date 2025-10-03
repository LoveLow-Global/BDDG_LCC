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
n           = 2^15
ps          = 0.01:0.01:0.5
n_replicate = 1000

results_file = "0_2to$(Int(log2(n)))_size_lcc.csv"

# Write header only once at the beginning
if !isfile(results_file)
    CSV.write(results_file, DataFrame(n=Int[], p=Float64[], rep=Int[], size_lcc=Int[]))
end

@time begin
    for p in ps
        # Create separate DataFrame for each thread
        results_per_thread = [DataFrame(n=Int[], p=Float64[], rep=Int[], size_lcc=Int[]) for _ in 1:Threads.nthreads()]
        
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
            size_lcc = nv(subg) # number of vertices of the largest SCC
            
            # Push the result to the DataFrame specific to the thread
            push!(results_per_thread[tid], (n, p, rep, size_lcc))
        end
        
        # 3. After the loop, merge the results from all threads.
        temp_results = vcat(results_per_thread...)
        
        # Write the combined results for this p to the file.
        CSV.write(results_file, temp_results; append=true)

        @info "Completed (n,p) = ($n, $p) at $(now())"
    end
end








