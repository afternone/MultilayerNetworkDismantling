## import packages
using LightGraphs
using GraphPlot
using DataStructures
using Plots
include("coreHLD.jl")
include("reverse_greedy.jl")

## compare HDA, HSDA and CoreHLD
layers = [erdos_renyi(10000, 1.5/10000) for _ in 1:2]
attack_nodes_HDA = HLDA(layers)
attack_nodes_HSDA = MHSDA(layers)
attack_nodes_mdecore = multilayer_decore(layers)
attack_nodes_treebreak, scomps = treebreak(layers, 100, attack_nodes_mdecore)

attack_nodes = vcat(attack_nodes_mdecore, attack_nodes_treebreak)
reinsert_nodes, presents = reverse_greedy!(layers, 100, attack_nodes)
g = decored_graph(layers, attack_nodes)
length.(connected_components(g)) |> maximum


y_HDA = recover_add_nodes(layers, attack_nodes_HDA)
y_HSDA = recover_add_nodes(layers, attack_nodes_HSDA)
y_mdecore = recover_add_nodes(layers, vcat(attack_nodes_mdecore, attack_nodes_treebreak))
nsum = sum(nv.(layers))
p = plot((1:length(y_HDA)) ./ nsum, y_HDA ./ nv(layers[1]), label="HLDA")
plot!(p,(1:length(y_HSDA)) ./ nsum, y_HSDA ./ nv(layers[1]), label="HSDA")
plot!(p,(1:length(y_mdecore)) ./ nsum, y_mdecore ./ nv(layers[1]), label="CoreHLD")
xlabel!(p, "fraction of removed nodes")
ylabel!(p, "fraction of nodes in GCC")
savefig(p, "result.png")

##
include("MultilayerNetworkDismantling.jl")

layers = [erdos_renyi(10000, 1.5/10000) for _ in 1:2]

attack_nodes_HLDA = MultilayerNetworkDismantling.HLDA(layers)
attack_nodes_CoreHLDA = MultilayerNetworkDismantling.CoreHLDA(layers)

y_HLDA = recover_add_nodes(layers, attack_nodes_HLDA)
y_CoreHLDA = recover_add_nodes(layers, attack_nodes_CoreHLDA)

nsum = sum(nv.(layers))
p = plot((1:length(y_HLDA)) ./ nsum, y_HLDA ./ nv(layers[1]), label="HLDA")
plot!(p, (1:length(y_CoreHLDA)) ./ nsum, y_CoreHLDA ./ nv(layers[1]), label="CoreHLDA")
xlabel!(p, "Fraction of nodes removed")
ylabel!(p, "Fraction of nodes in the LCC")
savefig(p, "result.png")
##