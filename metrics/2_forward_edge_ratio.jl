# Forward edge ratio computatoin doesn't take much time, plus this is not the main part of the research, so we do not optimize / multi-thread.
using Graphs, Graphs.SimpleGraphs
using SparseArrays, LinearAlgebra
using SimpleTraits, Statistics
using Base.Threads: @threads, nthreads, threadid
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

# PARAMETERS
n           = 2^10
ps          = 0.10:0.01:0.5
n_replicate = 10_000

results = DataFrame(n=Int[],p=Float64[], replicate=Int[], forward_edge_ratio=Float64[])

# Takes 1.32 to 1.76 seconds when n=2^15, ps=0.1,0.1,0.5, n_replicate=1
@time begin
    for p in ps
        for rep in 1:n_replicate
            g = ber_directed_divisor_graph(n, p)

            # Largest SCC
            scc_list   = strongly_connected_components_tarjan(g)
            comp_sizes = length.(scc_list)
            verts      = scc_list[argmax(comp_sizes)]
            subg, _    = induced_subgraph(g, verts)

            # compute forward-edge ratio
            fe = count(e -> src(e) < dst(e), edges(subg))
            tot = ne(subg)
            ratio = tot == 0 ? 0.0 : fe / tot

            push!(results, (n, p, rep, ratio))
        end
    end
end


if isfile("2_2to$(Int(log2(n)))_forward_edge_ratio.csv")
    CSV.write("2_2to$(Int(log2(n)))_forward_edge_ratio.csv", results; append=true)
else
    CSV.write("2_2to$(Int(log2(n)))_forward_edge_ratio.csv", results)
end