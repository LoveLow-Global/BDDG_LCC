using CSV
using DataFrames
using Statistics
using Plots
using Colors
using Compose

files = [
    "7_2to5_size_lcc_and_num_scc.csv",
    "7_2to10_size_lcc_and_num_scc.csv",
    "7_2to15_size_lcc_and_num_scc.csv",
    "7_2to20_size_lcc_and_num_scc.csv"
]

df = vcat([ CSV.read(f, DataFrame) for f in files ]...)

agg = combine(groupby(df, [:n, :p]), :size_plus_num_rate => mean => :mean_ratio)
sort!(agg, [:n, :p])

ns     = sort(unique(agg.n))
colors = palette(:tab10, length(ns))

# Plot
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
    # star at the Min (NOT Max)
    i_max = argmin(sub.mean_ratio)
    scatter!(plt,
             [sub.p[i_max]],
             [sub.mean_ratio[i_max]];
             marker     = :star5,
             markersize = 8,
             color      = col,
             label      = "")
end

xlabel!(plt, "p")
ylabel!(plt, "(Size of Largest SCC + Number of SCC -1)/n")
title!(plt, "Size of Largest SCC + Number of SCC -1")
display(plt)

savefig(plt, "7_size_lcc_and_num_scc_plot.pdf")