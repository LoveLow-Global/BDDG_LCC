using CSV
using DataFrames
using Statistics
using Plots

files = [
    "5_2to5_transitive_triangles.csv",
    "5_2to10_transitive_triangles.csv",
    "5_2to15_transitive_triangles.csv",
    "5_2to20_transitive_triangles.csv"
]

# read each into a DataFrame and vcat them together
df = vcat([ CSV.read(f, DataFrame) for f in files ]...)
# 1) Read & prepare
#df = CSV.read("5_transitive_triangles.csv", DataFrame)
df.r = df.tt ./ df.n
agg = combine(groupby(df, [:n, :p]), :r => mean => :mean_ratio)
sort!(agg, [:n, :p])

# 2) Build a palette & ordered list of n's
ns     = sort(unique(agg.n))
colors = palette(:tab10, length(ns))

# 3) Plot loop with pre-computed colors
plt = plot()
for (i, nval) in enumerate(ns)
    col = colors[i]
    sub = @view agg[agg.n .== nval, :]
    # line + circles
    plot!(plt,
          sub.p,
          sub.mean_ratio;
          label      = "n=2^$(Int(log2(nval)))",
          color      = col,
          marker     = :circle,
          markersize = 4,
          lw         = 2)
    # star at the max
    i_max = argmax(sub.mean_ratio)
    scatter!(plt,
             [sub.p[i_max]],
             [sub.mean_ratio[i_max]];
             marker     = :star5,
             markersize = 8,
             color      = col,
             label      = "")
end

# 4) Final tweaks
xlabel!(plt, "p")
ylabel!(plt, "Transitive Triangles / n")
title!(plt, "Transitive Triangles of the Largest SCC")

display(plt)

# savefig(plt, "5_Transative_triangles_plot.png")
using Compose
savefig(plt, "5_Transative_triangles_plot.pdf")