using LightGraphs
using Random
using DelimitedFiles
include("MultilayerNetworkDismantling.jl")

datafile = ARGS[1]
data = readdlm(datafile, Int)
n = maximum(data[:,2:3]) + 1
M = length(unique(data[:,1]))
layers = [SimpleGraph(n) for _ in 1:M]
for i in 1:size(data,1)
    l, u, v = data[i, 1:3]
    u != v && add_edge!(layers[l], u, v)
end
attacked_nodes = MultilayerNetworkDismantling.CoreHLDA(layers)
resultfile = splitext(datafile)[1]*".dismantling_set"
open(resultfile, "w") do f
    println(f, "layer", '\t', "id")
    for node in attacked_nodes
        println(f, node.layer, '\t', node.id)
    end
end