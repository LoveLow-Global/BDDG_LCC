using CSV
using DataFrames
using Statistics
using Plots
using Colors
using Compose

files = [
    "3_2to5_avg_path_length.csv",
    "3_2to10_avg_path_length.csv",
    "3_2to15_avg_path_length.csv",
    "3_2to20_avg_path_length.csv"
]

df = vcat([ CSV.read(f, DataFrame) for f in files ]...)

agg = combine(groupby(df, [:n, :p]), :apl => mean => :mean_ratio)

sort!(agg, [:n, :p])

ns     = sort(unique(agg.n)).
colors = palette(:tab10, length(ns))

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
    # Star at Max Not Required for Average Path Length, unless we are investigating in really small p.
    #=
    i_max = argmax(sub.mean_ratio)
    scatter!(plt,
             [sub.p[i_max]],
             [sub.mean_ratio[i_max]];
             marker     = :star5,
             markersize = 8,
             color      = col,
             label      = "")
    =#
end

xlabel!(plt, "p")
ylabel!(plt, "Avg Path Length")
title!(plt, "Avg Path length of the Largest SCC")
display(plt)

savefig(plt, "3_avg_path_length_plot.pdf")