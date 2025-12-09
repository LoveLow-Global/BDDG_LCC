using Plots
using Graphs
using GraphRecipes
using Random; Random.seed!(42)
using Compose

function ber_directed_divisibility_graph(n::Int, p::Float64)
    g = DiGraph(n)
    for i in 1:n, j in 2*i:i:n
        rand() < p ? add_edge!(g, i, j) : add_edge!(g, j, i)
    end
    return g
end

n, p = 8, 0.5
g    = ber_directed_divisibility_graph(n, p)

deg    = degree(g) # in+out degrees
labels = string.(1:nv(g))

plt = graphplot(
    g;
    node_weights = deg, # size ∝ degree
    nodesize     = 0.3, # maximum circle diameter
    method       = :spring,
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


savefig(plt, "Bernoulli-Directed Divisor Graph.pdf")
