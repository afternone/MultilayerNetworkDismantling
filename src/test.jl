## import packages
using Plots
using GraphPlot
include("MultilayerNetworkDismantling.jl")

layers = [erdos_renyi(500, 0.006) for _ in 1:2]
##
MultilayerNetworkDismantling.OAS(layers, 300)
##

attack_nodes_EMD = MultilayerNetworkDismantling.EMD(layers, num=300)
attack_nodes_HLD = MultilayerNetworkDismantling.HLD(layers)
attack_nodes_HLDA = MultilayerNetworkDismantling.HLDA(layers)
attack_nodes_HMD = MultilayerNetworkDismantling.HMD(layers)
attack_nodes_HMDA = MultilayerNetworkDismantling.HMDA(layers)
attack_nodes_HLCIA = MultilayerNetworkDismantling.HLCIA(layers, num=10000)
attack_nodes_HMCIA = MultilayerNetworkDismantling.HMCIA(layers, num=10000)
attack_nodes_CoreHLDA = MultilayerNetworkDismantling.CoreHLDA(layers)

##
y_EMD1 = Int[]
presents = [fill(true, nv(layers[1])) for _ in eachindex(layers)]
for (l,i) in attack_nodes_EMD
    g = MultilayerNetworkDismantling.merge_layer(layers, presents)
    push!(y_EMD1, MultilayerNetworkDismantling.LCC(g))
    presents[l][i] = false
end
##

y_EMD = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_EMD)
y_HLD = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLD)
y_HLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLDA)
y_HMD = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HMD)
y_HMDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HMDA)
y_HLCIA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLCIA)
y_HMCIA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HMCIA)
y_CoreHLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CoreHLDA)

nsum = sum(nv.(layers))
p = plot((1:length(y_HLD)) ./ nsum, y_HLD ./ nv(layers[1]), label="HLD")
plot!(p, (1:length(y_EMD)) ./ nsum, y_EMD ./ nv(layers[1]), label="EMD")
plot!(p, (1:length(y_HLDA)) ./ nsum, y_HLDA ./ nv(layers[1]), label="HLDA")
# plot!(p, (1:length(y_HMD)) ./ nsum, y_HMD ./ nv(layers[1]), label="HMD")
# plot!(p, (1:length(y_HMDA)) ./ nsum, y_HMDA ./ nv(layers[1]), label="HMDA")
# plot!(p, (1:length(y_HLCIA)) ./ nsum, y_HLCIA ./ nv(layers[1]), label="HLCIA")
# plot!(p, (1:length(y_HMCIA)) ./ nsum, y_HMCIA ./ nv(layers[1]), label="HMCIA")
plot!(p, (1:length(y_CoreHLDA)) ./ nsum, y_CoreHLDA ./ nv(layers[1]), label="CoreHLDA")
xlabel!(p, "Fraction of nodes removed")
ylabel!(p, "Fraction of nodes in the LCC")
#savefig(p, "result.png")
##
g1 = erdos_renyi(10000,1.5/10000)
g2 = erdos_renyi(10000,1.5/10000)
layers1 = MultilayerNetworkDismantling.correlated_multigraph(g1,g2,correlation="MP")
MultilayerNetworkDismantling.CoreHLDA(layers)