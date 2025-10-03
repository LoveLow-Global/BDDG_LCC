using Plots
using Graphs
using GraphRecipes      # brings in the graphplot recipe
using Random; Random.seed!(42)

# 1) build the Bernoulli‐directed divisibility graph
function ber_directed_divisibility_graph(n::Int, p::Float64)
    g = DiGraph(n)
    for i in 1:n, j in 2*i:i:n
        rand() < p ? add_edge!(g, i, j) : add_edge!(g, j, i)
    end
    return g
end

# 2) parameters
n, p = 8, 0.5
g    = ber_directed_divisibility_graph(n, p)

# 3) compute degrees, labels
deg    = degree(g)                   # in+out degrees
labels = string.(1:nv(g))            # ["1","2",…,"32"]

# 4) call graphplot with node_weights = deg
plt = graphplot(
    g;
    node_weights = deg,              # size ∝ degree
    nodesize     = 0.3,              # maximum circle diameter
    method       = :spring,          # force‐directed layout
    names        = labels,
    nodeshape    = :circle,
    nodecolor    = :turquoise3,
    edgecolor    = :black,
    curves       = true,
    arrow        = true,
    fontsize     = 10,
    title        = "Bernoulli-Directed Divisor Graph where G($n, $p)",
    background_color = :white
)

# 5) save as PDF
using Compose
savefig(plt, "Bernoulli-Directed Divisor Graph.pdf")
