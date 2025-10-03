using CSV
using DataFrames
using Statistics
using Plots
using Colors

# Read & prepare
files = [
    "4_2to5_diameter_results.csv",
    "4_2to10_diameter_results.csv",
    "4_2to15_diameter_results.csv",
    "4_2to20_diameter_results.csv"
]

# read each into a DataFrame and vcat them together
df = vcat([ CSV.read(f, DataFrame) for f in files ]...)
#df = CSV.read("4_diameter_results.csv", DataFrame)

agg = combine(groupby(df, [:n, :p]), :diameter => mean => :mean_ratio)

sort!(agg, [:n, :p])

# Build a palette & ordered list of n's
ns     = sort(unique(agg.n))
#colors = [ HSV(i/length(ns), 0.7, 0.9) for i in 1:length(ns) ]   # or palette(:tab10, length(ns)), etc.
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
ylabel!(plt, "Diameter")
title!(plt, "Diameter of the Largest SCC")

display(plt)

using Compose
savefig(plt, "4_diameter_plot.pdf")