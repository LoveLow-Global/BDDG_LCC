using Graphs, Graphs.SimpleGraphs
using SparseArrays, LinearAlgebra
using SimpleTraits, Statistics
using Base.Threads: @threads, nthreads, threadid
using Random;Random.seed!(42)
using Dates
using CSV, DataFrames

Threads.nthreads() # Check before running

function ber_directed_divisor_graph(n::Int, p::Float64; rng::AbstractRNG)
    g = DiGraph(n)
    @inbounds for i in 1:n, j in 2*i:i:n
        if rand(rng) < p
            add_edge!(g, i, j)
        else
            add_edge!(g, j, i)
        end 
    end
    return g
end

function diameter_exact(g::DiGraph)::Int
    n = nv(g)
    out_list = [collect(outneighbors(g, v)) for v in 1:n]
    in_list  = [collect(inneighbors(g, v))  for v in 1:n]

    active  = trues(n)
    visited = falses(n)
    queue   = Vector{Int}(undef, n)
    distbuf = fill(-1, n)
    diam    = 0

    vs = collect(1:n)
    sort!(vs, by = v -> -(length(out_list[v]) + length(in_list[v])))

    for u in vs
        if !active[u]; continue; end

        # Forward BFS
        fill!(visited, false)
        visited[u] = true
        queue[1]   = u
        front, back, level_end = 1, 2, 1
        e = 0
        while front < back
            v = queue[front]; front += 1
            @inbounds for w in out_list[v]
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

        # Backward BFS to prune
        dmax = diam - e
        if dmax ≥ 0
            fill!(distbuf, -1)
            distbuf[u] = 0
            queue[1]   = u
            front, back = 1, 2
            while front < back
                v = queue[front]; front += 1
                if distbuf[v] ≥ dmax
                    continue
                end
                @inbounds for w in in_list[v]
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

        if !any(active); break; end
    end
    return diam
end

# PARAMETERS
n           = 2^20
ps          = 0.01:0.01:0.50
n_replicate = 10

results_file = "4_2to$(Int(log2(n)))_diameter_results.csv"
if !isfile(results_file)
    CSV.write(results_file, DataFrame(n=Int[], p=Float64[], rep=Int[], diameter=Int[]))
end


@time begin
    for p in ps
        local_rows = Vector{NamedTuple{(:n,:p,:rep,:diameter),Tuple{Int,Float64,Int,Int}}}(undef, n_replicate)

        @threads for rep in 1:n_replicate
            g = ber_directed_divisibility_graph(n, p)

            # Largest SCC
            scc_list   = strongly_connected_components_tarjan(g)
            comp_sizes = length.(scc_list)
            verts      = scc_list[argmax(comp_sizes)]
            subg, _    = induced_subgraph(g, verts)

            d = diameter_exact(subg)
            local_rows[rep] = (n, p, rep, d)
        end

        CSV.write(results_file, DataFrame(local_rows); append=true)
        @info "Completed p = $p with $(n_replicate) reps @ $(now())"
    end
end