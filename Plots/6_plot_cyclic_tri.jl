using CSV
using DataFrames
using Statistics
using Plots

files = [
    "6_2to5_cyclic_triangles.csv",
    "6_2to10_cyclic_triangles.csv",
    "6_2to15_cyclic_triangles.csv",
    "6_2to20_cyclic_triangles.csv"
]

df = vcat([ CSV.read(f, DataFrame) for f in files ]...)

df.r = df.ct ./ df.n
agg = combine(groupby(df, [:n, :p]), :r => mean => :mean_ratio)
sort!(agg, [:n, :p])


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

        # Star at Max Not Required of Cyclic Triangles, it'll be 0.5 anyway.
end

xlabel!(plt, "p")
ylabel!(plt, "Cyclic Triangles / n")
title!(plt, "Cyclic Triangles of the Largest SCC")
display(plt)

using Compose
savefig(plt, "6_cyclic_triangles_plot.pdf")