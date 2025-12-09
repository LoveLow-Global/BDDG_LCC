using CSV
using DataFrames
using Statistics
using Plots
using Compose

files = [
    "0_2to5_size_lcc.csv",
    "0_2to10_size_lcc.csv",
    "0_2to15_size_lcc.csv",
    "0_2to20_size_lcc.csv"
]

df = vcat([ CSV.read(f, DataFrame) for f in files ]...)

df.r = df.size_lcc ./ df.n
agg = combine(groupby(df, [:n, :p]), :r => mean => :mean_ratio)
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

end

plot!(plt, ylim=(0, 1.0))
xlabel!(plt, "p")
ylabel!(plt, "Largest SCC size / n")
title!(plt, "Largest SCC per n")
display(plt)

savefig(plt, "0_size_largest_SCC_plot.pdf")