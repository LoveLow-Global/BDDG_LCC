using Graphs, Graphs.SimpleGraphs
using SparseArrays, LinearAlgebra
using SimpleTraits, Statistics
using Base.Threads: @threads, nthreads, threadid
using Random; Random.seed!(42)
using Dates
using CSV, DataFrames

# Check
Threads.nthreads()

# Divisibility graph with Bernoulli orientation (global RNG)
function ber_directed_divisor_graph(n::Int, p::Float64)
    g = DiGraph(n)
    @inbounds for i in 1:n, j in 2*i:i:n
        if rand() < p
            add_edge!(g, i, j)
        else
            add_edge!(g, j, i)
        end
    end
    return g
end

# Build CSR for induced subgraph on verts (local vertices are 1..n2)
# Returns (rowptr_out::Vector{Int}, col_out::Vector{Int}, rowptr_in::Vector{Int}, col_in::Vector{Int})
function build_csr_induced(g::DiGraph, verts::Vector{Int})
    n2 = length(verts)
    # global->local map (0 means absent)
    gl2lc = zeros(Int, nv(g))
    @inbounds for (i, v) in enumerate(verts)
        gl2lc[v] = i
    end

    # Pass 1: count out-degrees within the set
    outdeg = zeros(Int, n2)
    @inbounds for (i, v) in enumerate(verts)
        cnt = 0
        for w in outneighbors(g, v)
            if gl2lc[w] != 0
                cnt += 1
            end
        end
        outdeg[i] = cnt
    end

    # Prefix sums -> rowptr_out (1-based CSR)
    rowptr_out = Vector{Int}(undef, n2 + 1)
    rowptr_out[1] = 1
    @inbounds for i in 1:n2
        rowptr_out[i+1] = rowptr_out[i] + outdeg[i]
    end
    col_out = Vector{Int}(undef, rowptr_out[end] - 1)

    # Pass 2: fill col_out and accumulate indegrees
    pos_out = copy(rowptr_out)
    indeg   = zeros(Int, n2)
    @inbounds for (i, v) in enumerate(verts)
        p = pos_out[i]
        for w in outneighbors(g, v)
            j = gl2lc[w]
            if j != 0
                col_out[p] = j
                p += 1
                indeg[j] += 1
            end
        end
        pos_out[i] = p
    end

    # Build CSR for incoming edges from counts
    rowptr_in = Vector{Int}(undef, n2 + 1)
    rowptr_in[1] = 1
    @inbounds for i in 1:n2
        rowptr_in[i+1] = rowptr_in[i] + indeg[i]
    end
    col_in = Vector{Int}(undef, rowptr_in[end] - 1)

    # Fill incoming CSR by traversing out CSR once
    pos_in = copy(rowptr_in)
    @inbounds for i in 1:n2
        s = rowptr_out[i]
        e = rowptr_out[i+1] - 1
        for k in s:e
            j = col_out[k]
            col_in[pos_in[j]] = i
            pos_in[j] += 1
        end
    end

    return rowptr_out, col_out, rowptr_in, col_in
end

# Helpers to get degree from CSR
@inline function outdeg_csr(rowptr::Vector{Int}, v::Int)
    return rowptr[v+1] - rowptr[v]
end
@inline function indeg_csr(rowptr::Vector{Int}, v::Int)
    return rowptr[v+1] - rowptr[v]
end

# Exact diameter on CSR adjacency with epoch trick (no repeated clears)
function diameter_exact_csr(rpo::Vector{Int}, co::Vector{Int},
                            rpi::Vector{Int}, ci::Vector{Int})::Int
    n = length(rpo) - 1
    active    = fill(true, n)              # Vector{Bool}
    remaining = n

    visited_epoch = zeros(Int, n)          # valid when == token
    dist_epoch    = zeros(Int, n)          # valid when == dtok
    dist          = fill(-1, n)
    queue         = Vector{Int}(undef, n)

    token = 0
    dtok  = 0
    diam  = 0

    # High-(in+out)-degree first
    degsum = Vector{Int}(undef, n)
    @inbounds for v in 1:n
        degsum[v] = outdeg_csr(rpo, v) + indeg_csr(rpi, v)
    end
    vs = collect(1:n)
    sort!(vs, by = v -> -degsum[v])

    @inbounds for u in vs
        if !active[u]; continue; end

        # Forward BFS for ecc(u)
        token += 1
        visited_epoch[u] = token
        queue[1] = u
        front = 1; back = 2; level_end = 1
        e = 0
        while front < back
            v = queue[front]; front += 1
            s = rpo[v]; eoi = rpo[v+1] - 1
            for k in s:eoi
                w = co[k]
                if visited_epoch[w] != token
                    visited_epoch[w] = token
                    queue[back] = w; back += 1
                end
            end
            if front > level_end && front < back
                e += 1
                level_end = back - 1
            end
        end
        if e > diam; diam = e; end

        # Depth-limited backward BFS to prune
        dmax = diam - e
        if dmax >= 0
            dtok += 1
            dist[u] = 0; dist_epoch[u] = dtok
            queue[1] = u
            front = 1; back = 2
            while front < back
                v = queue[front]; front += 1
                if dist[v] >= dmax; continue; end
                s = rpi[v]; eii = rpi[v+1] - 1
                for k in s:eii
                    w = ci[k]
                    if dist_epoch[w] != dtok
                        dist_epoch[w] = dtok
                        dist[w] = dist[v] + 1
                        queue[back] = w; back += 1
                    end
                end
            end
            # prune
            for v in 1:n
                if active[v] && dist_epoch[v] == dtok && (dist[v] + e) <= diam
                    active[v] = false
                    remaining -= 1
                end
            end
        end

        if remaining == 0; break; end
    end
    return diam
end

n           = 2^20
ps          = 0.5:-0.01:0.45
n_replicate = 10

results_file = "4_2to$(Int(log2(n)))_diameter_results_ALT_ALT.csv"
if !isfile(results_file)
    CSV.write(results_file, DataFrame(n=Int[], p=Float64[], rep=Int[], diameter=Int[]))
end

@time begin
    for p in ps
        local_rows = Vector{NamedTuple{(:n,:p,:rep,:diameter),Tuple{Int,Float64,Int,Int}}}(undef, n_replicate)

        @threads for rep in 1:n_replicate
            g = ber_directed_divisor_graph(n, p)

            # L-SCC
            scc_list   = strongly_connected_components_tarjan(g)
            comp_sizes = length.(scc_list)
            verts      = scc_list[argmax(comp_sizes)]

            # CSR adjacencies (no push!, contiguous)
            rpo, co, rpi, ci = build_csr_induced(g, verts)

            # exact diameter
            d = diameter_exact_csr(rpo, co, rpi, ci)

            local_rows[rep] = (n, p, rep, d)
        end

        CSV.write(results_file, DataFrame(local_rows); append=true)
        @info "Completed p = $p with $(n_replicate) reps @ $(now())"
    end
end
 