# Forward edge ratio computatoin doesn't take much time, so we do not optimize / multi-thread.
using Graphs, Graphs.SimpleGraphs
using SparseArrays, LinearAlgebra
using SimpleTraits, Statistics
using Base.Threads: @threads, nthreads, threadid, Atomic, atomic_add!
using Random; Random.seed!(42)
using Dates
using CSV, DataFrames

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

function bfs_distances(g::DiGraph, src::Int)
    n = nv(g)
    dist = fill(-1, n)
    dist[src] = 0
    q = Vector{Int}(undef, n)
    front, back = 1, 1
    q[back] = src; back += 1

    while front < back
        u = q[front]; front += 1
        for v in outneighbors(g, u)
            if dist[v] < 0
                dist[v] = dist[u] + 1
                q[back] = v; back += 1
            end
        end
    end

    return dist
end

function avg_path_length(subg::DiGraph; samples::Int)
    n_sub = nv(subg)
    dists = Float64[]
    for _ in 1:samples
        u, v = rand(1:n_sub), rand(1:n_sub)
        dist = bfs_distances(subg, u) # dist is a Vector
        push!(dists, dist[v])         # index it directly
    end
    return isempty(dists) ? 0.0 : mean(dists)
end

# PARAMETERS
n           = 2^20
ps          = 0.01:0.01:0.5
n_replicate = 5

results_file = "3_2to$(Int(log2(n)))_avg_path_length.csv"

if !isfile(results_file)
    CSV.write(results_file, DataFrame(n=Int[], p=Float64[], rep=Int[], apl=Float64[]))
end


# 14 sec when n=2^15, ps=0.1,0.1,0.5, n_replicate=1
@time begin
    for p in ps
        temp_results = DataFrame(n=Int[],p=Float64[], rep=Int[], apl=Float64[])
        
        for rep in 1:n_replicate
            g = ber_directed_divisor_graph(n, p)

            # Largest SCC
            scc_list   = strongly_connected_components_tarjan(g)
            comp_sizes = length.(scc_list)
            verts      = scc_list[argmax(comp_sizes)]
            subg, _    = induced_subgraph(g, verts)

            # average path-length with 1000 samples
            apl = avg_path_length(subg; samples=1000) # Edit sample size if needed

            push!(temp_results, (n, p, rep, apl))
        end
        
        CSV.write(results_file, temp_results; append=true)

        @info "Completed (n,p) = ($n, $p) at $(now())"
    end
end