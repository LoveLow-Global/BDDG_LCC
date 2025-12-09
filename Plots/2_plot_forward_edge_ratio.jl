using CSV
using DataFrames
using Statistics
using Plots
using Compose

files = [
    "2_2to5_forward_edge_ratio.csv",
    "2_2to10_forward_edge_ratio.csv",
    "2_2to15_forward_edge_ratio.csv",
    "2_2to20_forward_edge_ratio.csv"
]

df = vcat([ CSV.read(f, DataFrame) for f in files ]...)

agg = combine(groupby(df, [:n, :p]),
              :forward_edge_ratio => mean => :mean_ratio)

sort!(agg, [:n, :p])

plt = plot()

for nval in unique(agg.n)
    sub = @view agg[agg.n .== nval, :]
    plot!(plt,
          sub.p,
          sub.mean_ratio;
          label="n=2^$(Int(log2(nval)))",
          marker=:circle,
          linestyle=:auto)
end

xlabel!("p")
ylabel!("Forward Edge Ratio")
title!("Forward Edge Ratio vs. p for Different n")
display(plt)

savefig(plt, "2_forward_edge_ratio.pdf")