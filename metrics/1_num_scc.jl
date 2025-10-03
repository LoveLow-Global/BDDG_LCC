# num_scc.jl

using Graphs, Graphs.SimpleGraphs
using SparseArrays, LinearAlgebra
using SimpleTraits, Statistics
using Base.Threads: @threads, nthreads, threadid
using Random; Random.seed!(42)
using Dates
using CSV, DataFrames

Threads.nthreads()

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
ps          = 0.01:0.01:0.5
n_replicate = 10_000

results = DataFrame(n=Int[],p=Float64[], replicate=Int[], num_scc=Int[])

# 0.56 to 0.95 sec when n=2^15, ps=0.1,0.1,0.5, n_replicate=1
@time begin
    for p in ps
    for rep in 1:n_replicate
        g = ber_directed_divisor_graph(n, p)
        scc_list = strongly_connected_components_tarjan(g)
        push!(results, (n, p, rep, length(scc_list)))
    end
end
end


if isfile("1_2to$(Int(log2(n)))_num_scc_results.csv")
    CSV.write("1_2to$(Int(log2(n)))_num_scc_results.csv", results; append=true)
else
    CSV.write("1_2to$(Int(log2(n)))_num_scc_results.csv", results)  # write with header
end