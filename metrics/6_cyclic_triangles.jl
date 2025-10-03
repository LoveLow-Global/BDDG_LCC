# cyclic_triangles.jl

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

# simple 3‐cycle count
#=
function directed_cyclic_triangles(g::DiGraph)::Int
    A = adjacency_matrix(g)
    return tr(A * A * A) ÷ 6
end
=#

function directed_cyclic_triangles(g::DiGraph)::Int
    n = nv(g)
    mark = falses(n)          # one Bool per vertex for marking
    total = 0

    for u in vertices(g)
        # 1) mark all inneighbors of u (those w with w → u)
        for w in inneighbors(g, u)
            mark[w] = true
        end

        # 2) for each two‐step path u → v → w, check if we have w → u
        for v in outneighbors(g, u)
            for w in outneighbors(g, v)
                @inbounds if mark[w]
                    total += 1
                end
            end
        end

        # 3) unmark
        for w in inneighbors(g, u)
            mark[w] = false
        end
    end

    # each directed 3-cycle shows up 3 times in these walks,
    # since u → v → w → u, v → w → u → v, w → u → v → w
    # so divide by 3 to get the unique count
    return total ÷ 3
end

# PARAMETERS
n           = 2^10
ps          = 0.01:0.01:0.5
n_replicate = 10000

results = DataFrame(n=Int[],p=Float64[], replicate=Int[], cyclic_triangles=Int[])

results_file = "6_2to$(Int(log2(n)))_cyclic_triangles.csv"

if !isfile(results_file)
    CSV.write(results_file, DataFrame(n=Int[], p=Float64[], rep=Int[], ct=Float64[]))
end

#
@time begin
    for p in ps
        for rep in 1:n_replicate
            g = ber_directed_divisor_graph(n, p)
            # find largest SCC
            scc_list   = strongly_connected_components_tarjan(g)
            comp_sizes = length.(scc_list)
            verts      = scc_list[argmax(comp_sizes)]
            subg, _    = induced_subgraph(g, verts)

            # compute diameter
            ct = sum(directed_cyclic_triangles(subg))

            # push n, p, rep, and diameter
            push!(results, (n, p, rep, ct))
        end

        # Write results for this p
        CSV.write(results_file, results; append=true)

        @info "Completed (n,p) = ($n, $p) at $(now())"
    end
end
