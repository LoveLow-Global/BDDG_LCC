using CSV
using DataFrames
using Statistics
using Plots
using Compose

files = [
    "1_2to5_num_scc_results.csv",
    "1_2to10_num_scc_results.csv",
    "1_2to15_num_scc_results.csv",
    "1_2to20_num_scc_results.csv"
]

df = vcat([ CSV.read(f, DataFrame) for f in files ]...)
df.r = df.num_scc ./ df.n
agg = combine(groupby(df, [:n, :p]), :r => mean => :mean_ratio)
sort!(agg, [:n, :p])

ns     = sort(unique(agg.n))
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

end

xlabel!(plt, "p")
ylabel!(plt, "Number of SCCs / n")
title!(plt, "Number of SCCs per n")
display(plt)

savefig(plt, "1_plot_num_scc.pdf")