# diameter.jl
using Graphs, Graphs.SimpleGraphs
using SparseArrays, LinearAlgebra
using SimpleTraits, Statistics
using Base.Threads: @threads, nthreads, threadid
using Random; Random.seed!(42)
using Dates
using CSV, DataFrames

function ber_directed_divisibility_graph(n::Int64, p::Float64)
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

#=
# the same double‐sweep BFS diameter as in metrics_with_opt_diam.jl
function approximate_diameter(g::DiGraph; num_samples::Int=10000)
    n = nv(g)
    if n == 0
        return 0
    end

    # pick a few random start‐vertices (plus 1 and n)
    starts = Set{Int}()
    if n >= 1
        push!(starts, 1, n)
    end
    while length(starts) < num_samples + 2 && length(starts) < n
        push!(starts, rand(1:n))
    end

    # inner BFS that returns (distance, farthest_vertex)
    function bfs_farthest(src::Int)
        visited = falses(n)
        visited[src] = true
        queue = Vector{Int}(undef, n)
        front, back = 1, 1
        queue[back] = src; back += 1

        dist = 0
        level_end = 1
        farthest_vertex = src

        while front < back
            v = queue[front]; front += 1
            for w in outneighbors(g, v)
                if !visited[w]
                    visited[w] = true
                    queue[back] = w; back += 1
                end
            end
            if front > level_end
                dist += 1
                level_end = back - 1
                farthest_vertex = queue[level_end]
            end
        end

        return dist, farthest_vertex
    end

    # 1st sweep over all start vertices
    best_distance, best_vertex = -1, 1
    for s in starts
        d, v_far = bfs_farthest(s)
        if d > best_distance
            best_distance, best_vertex = d, v_far
        end
    end

    # 2nd sweep from the farthest found
    d2, _ = bfs_farthest(best_vertex)
    best_distance = max(best_distance, d2)

    return best_distance
end
=#

function diameter_exact2(g::DiGraph)::Int
    n = nv(g)
    # cache neighbors
    out_list = [collect(outneighbors(g, v)) for v in 1:n]
    in_list  = [collect(inneighbors(g, v))  for v in 1:n]

    active  = trues(n)
    visited = falses(n)
    queue   = Vector{Int}(undef, n)
    distbuf = fill(-1, n)
    diam    = 0

    # process vertices in descending total‐degree order
    vs = collect(1:n)
    sort!(vs, by = v -> -(length(out_list[v]) + length(in_list[v])))

    for u in vs
        if !active[u]
            continue
        end

        # forward BFS from u → compute e = ecc(u)
        fill!(visited, false)
        visited[u] = true
        queue[1]   = u
        front, back, level_end = 1, 2, 1
        e = 0
        while front < back
            v = queue[front]; front += 1
            for w in out_list[v]
                if !visited[w]
                    visited[w] = true
                    queue[back] = w
                    back += 1
                end
            end
            if front > level_end && front < back
                e += 1
                level_end = back - 1
            end
        end
        diam = max(diam, e)

        # depth‐limited backward BFS from u → prune
        dmax = diam - e
        if dmax ≥ 0
            fill!(distbuf, -1)
            distbuf[u] = 0
            queue[1]    = u
            front, back = 1, 2
            while front < back
                v = queue[front]; front += 1
                if distbuf[v] ≥ dmax
                    continue
                end
                for w in in_list[v]
                    if distbuf[w] < 0
                        distbuf[w] = distbuf[v] + 1
                        queue[back] = w
                        back += 1
                    end
                end
            end
            @inbounds for v in 1:n
                if active[v] && distbuf[v] ≥ 0 && distbuf[v] + e ≤ diam
                    active[v] = false
                end
            end
        end

        if !any(active)
            break
        end
    end

    return diam
end

# PARAMETERS
n           = 2^15
ps          = 0.01:0.01:0.50
n_replicate = 100

results = DataFrame(n=Int[],p=Float64[], replicate=Int[], diameter=Int[])

results_file = "4_2to$(Int(log2(n)))_diameter_results_ALT.csv"

if !isfile(results_file)
    CSV.write(results_file, DataFrame(n=Int[], p=Float64[], rep=Int[], diameter=Float64[]))
end

# 213 sec when n=2^15, ps=0.1:0.1:0.5, n_replicate=1
# 1739 sec when n=2^20, ps=0.1, n_replicate=1
@time begin
    for p in ps
        # collect only this p’s results
        local_results = DataFrame(n=Int[], p=Float64[], rep=Int[], diameter=Int[])
        for rep in 1:n_replicate
            g = ber_directed_divisibility_graph(n, p)
            # find largest SCC
            scc_list   = strongly_connected_components_tarjan(g)
            comp_sizes = length.(scc_list)
            verts      = scc_list[argmax(comp_sizes)]
            subg, _    = induced_subgraph(g, verts)

            # compute diameter
            d = diameter_exact2(subg)

            # push only into this p’s DataFrame
            push!(local_results, (n, p, rep, d))
        end

        # append *only* this p’s block to file
        CSV.write(results_file, local_results; append=true)
        @info "Completed p = $p at $(now())"
    end
end








