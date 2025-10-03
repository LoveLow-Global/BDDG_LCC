using CSV
using DataFrames
using Statistics
using Plots



# 1) Read & prepare
files = [
    "1_2to5_num_scc_results.csv",
    "1_2to10_num_scc_results.csv",
    "1_2to15_num_scc_results.csv",
    "1_2to20_num_scc_results.csv"
]

# read each into a DataFrame and vcat them together
df = vcat([ CSV.read(f, DataFrame) for f in files ]...)
#df = CSV.read("1_num_scc_results.csv", DataFrame)
df.r = df.num_scc ./ df.n
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
ylabel!(plt, "Number of SCCs / n")
title!(plt, "Number of SCCs per n")

display(plt)

# savefig(plt, "6_cyclic_triangles_plot.png")
using Compose
savefig(plt, "1_plot_num_scc.pdf")