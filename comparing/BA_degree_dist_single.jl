# Number of SCCs and $N$



using Graphs
#using GraphRecipes
using Plots
using StatsBase
using Statistics

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

# just input (n,p) here
g = ber_directed_divisor_graph(2^20, 0.5)

# find largest SCC
scc_list   = strongly_connected_components_tarjan(g)
comp_sizes = length.(scc_list)
verts      = scc_list[argmax(comp_sizes)]
subg, _    = induced_subgraph(g, verts)

vertix_num = nv(subg)
edge_num = ne(subg)

avg_degree = edge_num / vertix_num

BA_graph = barabasi_albert(vertix_num, Int(round(avg_degree)))


# Calculate degree distribution for BDDG
BDDG_degrees = degree(subg)
BDDG_degree_hist = countmap(BDDG_degrees)
BDDG_degree_values = collect(keys(BDDG_degree_hist))  # Different degree values (k)
BDDG_frequencies = collect(values(BDDG_degree_hist))  # Number of nodes with each degree

# Calculate degree distribution for BA
BA_degrees = degree(BA_graph)
BA_degree_hist = countmap(BA_degrees)
BA_degree_values = collect(keys(BA_degree_hist))  # Different degree values (k)
BA_frequencies = collect(values(BA_degree_hist))  # Number of nodes with each degree


# Plot degree dist of Bernoulli-Directed Divisor Graph
plt = scatter(BDDG_degree_values, BDDG_frequencies,
    xlab="Number of Edges (k)", 
    ylab="Number of Vertices with k Edges", 
    title=#"Degree Distribution Comparison - BDDG vs BA",
    label=#"BDDG Model Degree Distribution",
    legend=:topright,
    α= 0.5,
    markersize=2,
    markercolor=:black, # slateblue1
    markerstrokecolor=:slateblue1,
    xscale=:log10,
    yscale=:log10)

# Plot degree dist of BA model graph
scatter!(plt,BA_degree_values, BA_frequencies,
    #label="BA Model Degree Distribution",
    α = 0.5,
    markersize=2,
    markercolor=:yellow,
    markerstrokecolor=:yellow)

using Compose
savefig(plt, "BA_degree_dist_plot.pdf")