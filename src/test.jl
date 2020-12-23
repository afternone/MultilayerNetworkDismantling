## import packages
using Plots
using GraphPlot
using LightGraphs
include("MultilayerNetworkDismantling.jl")
##
n = 10000 # number of nodes
k = 3 # average degree
M = 2 # number of layers
layers = [erdos_renyi(n, k/n) for _ in 1:M]
##
#MultilayerNetworkDismantling.OAS(layers, 300)
#attack_nodes_HLD = MultilayerNetworkDismantling.HLD(layers)
attack_nodes_EMD = MultilayerNetworkDismantling.EMD(layers, num=div(nv(layers[1]),2))
attack_nodes_HLDA = MultilayerNetworkDismantling.HLDA(layers)
#attack_nodes_HMD = MultilayerNetworkDismantling.HMD(layers)
attack_nodes_HMDA = MultilayerNetworkDismantling.HMDA(layers)
attack_nodes_CI = MultilayerNetworkDismantling.HLCIA(layers, num=10000)
#attack_nodes_HMCIA = MultilayerNetworkDismantling.HMCIA(layers, num=10000)
attack_nodes_CoreHLDA = MultilayerNetworkDismantling.CoreHLDA(layers, 2)
##
y_EMD = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_EMD)
#y_HLD = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLD)
y_HLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLDA)
#y_HMD = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HMD)
y_HMDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HMDA)
y_CI = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CI)
#y_HMCIA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HMCIA)
y_CoreHLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CoreHLDA)
## plot
nsum = sum(nv.(layers))
pyplot()
#p_HLD = plot((1:length(y_HLD)) ./ nsum, y_HLD ./ nv(layers[1]), label="HLD")
p_EMD = plot((1:length(y_EMD)) ./ nsum, y_EMD ./ nv(layers[1]), label="EMD")
p_HLDA = plot((1:length(y_HLDA)) ./ nsum, y_HLDA ./ nv(layers[1]), label="HLDA")
# plot!(p, (1:length(y_HMD)) ./ nsum, y_HMD ./ nv(layers[1]), label="HMD")
p_HMDA = plot((1:length(y_HMDA)) ./ nsum, y_HMDA ./ nv(layers[1]), label="HMDA")
p_CI = plot((1:length(y_CI)) ./ nsum, y_CI ./ nv(layers[1]), label="CI")
# plot!(p, (1:length(y_HMCIA)) ./ nsum, y_HMCIA ./ nv(layers[1]), label="HMCIA")
p_CoreHLDA = plot((1:length(y_CoreHLDA)) ./ nsum, y_CoreHLDA ./ nv(layers[1]), label="CoreHLDA")
p = plot(p_EMD, p_HLDA, p_HMDA, p_CI, p_CoreHLDA)
xlabel!(p, "Fraction of removed nodes")
ylabel!(p, "Relative size of the largest cluster")
#savefig(p, "result.png")
## correlated multiplex network
g1 = erdos_renyi(10000,1.5/10000)
g2 = erdos_renyi(10000,1.5/10000)
layers1 = MultilayerNetworkDismantling.correlated_multigraph(g1,g2,correlation="MP")
MultilayerNetworkDismantling.CoreHLDA(layers)
##