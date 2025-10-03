# transitive_triangles.jl

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
ps          = 0.01:0.01:0.50
n_replicate = 10000

results = DataFrame(n=Int[], p=Float64[], replicate=Int[], transitive_triangles=Int[])

results_file = "5_2to$(Int(log2(n)))_transitive_triangles.csv"

if !isfile(results_file)
    CSV.write(results_file, DataFrame(n=Int[], p=Float64[], rep=Int[], tt=Float64[]))
end

# 2.85 to 3 sec when n=2^15, ps=0.1,0.1,0.5, n_replicate=1
#=
@time begin
    for p in ps
        for rep in 1:n_replicate
            g = ber_directed_divisibility_graph(n, p)
            # find largest SCC
            scc_list   = strongly_connected_components_tarjan(g)
            comp_sizes = length.(scc_list)
            verts      = scc_list[argmax(comp_sizes)]
            subg, _    = induced_subgraph(g, verts)

            # compute transitive triangles
            tt = sum(triangles(subg))

            push!(results, (n, p, rep, tt))
        end
        # Write results for this p
        CSV.write(results_file, results; append=true)

        @info "Completed (n,p) = ($n, $p) at $(now())"
    end
end
=#


@time begin
    for p in ps
        for rep in 1:n_replicate
            g = ber_directed_divisor_graph(n, p)
            # find largest SCC
            scc_list   = strongly_connected_components_tarjan(g)
            comp_sizes = length.(scc_list)
            verts      = scc_list[argmax(comp_sizes)]
            subg, _    = induced_subgraph(g, verts)

            # compute transitive triangles
            tt = sum(triangles(subg))

            push!(results, (n, p, rep, tt))
        end
        let start_idx = nrow(results) - n_replicate + 1
            CSV.write(results_file, results[start_idx:end, :]; append=true)
        end
        @info "Completed (n,p) = ($n, $p) at $(now())"
    end
end

#CSV.write("transitive_triangles_results.csv", results)
#=
if isfile("5_transitive_triangles.csv")
    CSV.write("5_transitive_triangles.csv", df; append=true)
else
    CSV.write("5_transitive_triangles.csv", df)  # write with header
end
=#