## import packages
using LightGraphs
using GraphPlot
using DataStructures
using Plots
include("coreHLD.jl")

## compare HDA, HSDA and CoreHLD
layers = [erdos_renyi(10000, 1.5 / 10000) for _ in 1:2]
attack_nodes_HDA = MHDA(layers)
attack_nodes_HSDA = MHSDA(layers)
attack_nodes_mdecore = multilayer_decore(layers)
attack_nodes_treebreak, scomps = treebreak(layers, 1, attack_nodes_mdecore)

y_HDA = recover_add_nodes(layers, attack_nodes_HDA)
y_HSDA = recover_add_nodes(layers, attack_nodes_HSDA)
y_mdecore = recover_add_nodes(layers, vcat(attack_nodes_mdecore, attack_nodes_treebreak))
nsum = sum(nv.(layers))
p = plot((1:length(y_HDA)) ./ nsum, y_HDA ./ nv(layers[1]), label="HDA")
plot!(p,(1:length(y_HSDA)) ./ nsum, y_HSDA ./ nv(layers[1]), label="HSDA")
plot!(p,(1:length(y_mdecore)) ./ nsum, y_mdecore ./ nv(layers[1]), label="Mdecore")
xlabel!(p, "fraction of removed nodes")
ylabel!(p, "fraction of nodes in GCC")
savefig(p, "result.png")

##