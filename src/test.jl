## import packages
using Plots
using Plots.PlotMeasures
using GraphPlot
using LightGraphs
using JLD2
using Random
using LaTeXStrings
using DelimitedFiles
include("MultilayerNetworkDismantling.jl")
##
n = 10000
k = 3
g1 = erdos_renyi(n,k/n)
g2 = erdos_renyi(n,k/n)
layers = MultilayerNetworkDismantling.correlated_multigraph(g1,g2,correlation="UC")
##
g1 = barabasi_albert(n,div(k,2))
g2 = barabasi_albert(n,div(k,2))
layers = MultilayerNetworkDismantling.correlated_multigraph(g1,g2,correlation="UC")
##
attack_nodes_HLD = MultilayerNetworkDismantling.HLD(layers)
attack_nodes_HLDA = MultilayerNetworkDismantling.HLDA(layers)
attack_nodes_HALDA = MultilayerNetworkDismantling.HALDA(layers)
attack_nodes_HLCIA = MultilayerNetworkDismantling.HLCIA(layers)
attack_nodes_HACILDA = MultilayerNetworkDismantling.HACILDA(layers)
attack_nodes_CoreHLDA = MultilayerNetworkDismantling.CoreHLDA(layers)
##
GCC_HLD = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLD)
GCC_HLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLDA)
GCC_HALDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HALDA)
GCC_HLCIA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLCIA)
GCC_HACILDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HACILDA)
GCC_CoreHLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CoreHLDA)
##
i_HLD = findfirst(x->x≤100, GCC_HLD)
i_HLDA = findfirst(x->x≤100, GCC_HLDA)
i_HALDA = findfirst(x->x≤100, GCC_HALDA)
i_HLCIA = findfirst(x->x≤100, GCC_HLCIA)
i_HACILDA = findfirst(x->x≤100, GCC_HACILDA)
y_HLD = GCC_HLD[1:i_HLD]
y_HLDA = GCC_HLDA[1:i_HLDA]
y_HALDA = GCC_HALDA[1:i_HALDA]
y_HLCIA = GCC_HLCIA[1:i_HLCIA]
y_HACILDA = GCC_HACILDA[1:i_HACILDA]
y_CoreHLDA = GCC_CoreHLDA

p = plot((1:length(y_HLDA))./20000, y_HLDA./10000, label="HLDA")
plot!(p, (1:length(y_HLD))./20000, y_HLD./10000, label="HLD")
plot!(p, (1:length(y_HALDA))./20000, y_HALDA./10000, label="HALDA")
plot!(p, (1:length(y_HLCIA))./20000, y_HLCIA./10000, label="HLCIA")
plot!(p, (1:length(y_HACILDA))./20000, y_HACILDA./10000, label="HACILDA")
plot!(p, (1:length(y_CoreHLDA))./20000, y_CoreHLDA./10000, label="CoreHLDA")
plot!(p, legend=:bottomleft)
##
n = 10000
degs = 3:12
fc_CI = zeros(10,10)
fc_EMD = zeros(10)
fc_HLDA = zeros(10,10)
fc_HMDA = zeros(10,10)
fc_CoreHLDA = zeros(10,10)
for i in 1:10
    println(i)
    for j in eachindex(degs)
        g1 = erdos_renyi(n, degs[j]/n)
        g2 = erdos_renyi(n, degs[j]/n)
        layers = MultilayerNetworkDismantling.correlated_multigraph(g1,g2,correlation="UC")

        # attack_nodes_CI = MultilayerNetworkDismantling.HLCIA(layers, num=15000)
        # attack_nodes_EMD = MultilayerNetworkDismantling.EMD(layers, num=8000)
        # attack_nodes_HLDA = MultilayerNetworkDismantling.HLDA(layers)
        # attack_nodes_HMDA = MultilayerNetworkDismantling.HMDA(layers)
        attack_nodes_CoreHLDA = MultilayerNetworkDismantling.CoreHLDA(layers)

        # GCC_CI = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CI)
        # GCC_EMD = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_EMD)
        # GCC_HLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLDA)
        # GCC_HMDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HMDA)
        GCC_CoreHLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CoreHLDA)

        # i_CI = findfirst(x->x≤100, GCC_CI)
        # i_EMD = findfirst(x->x≤100, GCC_EMD)
        # i_HLDA = findfirst(x->x≤100, GCC_HLDA)
        # i_HMDA = findfirst(x->x≤100, GCC_HMDA)
        i_CoreHLDA = length(GCC_CoreHLDA)
        
        nsum = sum(nv.(layers))
        # fc_CI[i,j] = i_CI/nsum
        # fc_EMD[j] = i_EMD/nsum
        # fc_HLDA[i,j] = i_HLDA/nsum
        # fc_HMDA[i,j] = i_HMDA/nsum
        fc_CoreHLDA[i,j] = i_CoreHLDA/nsum
    end
end
@save "CoreHLDA_ER_UC_ER_k_3_12.jld2" fc_CoreHLDA
## time complexity analysis
n = [2^i for i in 10:19]
T = zeros(10,10)
fc_CoreHLDA = zeros(10,10)
for i in 1:10
    for j in eachindex(n)
        println(i,'\t', j)
        g1 = erdos_renyi(n[j], 6/n[j])
        g2 = erdos_renyi(n[j], 6/n[j])
        layers = MultilayerNetworkDismantling.correlated_multigraph(g1,g2,correlation="UC")
        T[i,j] = @elapsed attack_nodes_CoreHLDA = MultilayerNetworkDismantling.CoreHLDA(layers)
        GCC_CoreHLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CoreHLDA)
        i_CoreHLDA = length(GCC_CoreHLDA)
        nsum = sum(nv.(layers))
        fc_CoreHLDA[i,j] = i_CoreHLDA/nsum
    end
end
## plot the runtime
gr(size=(400,300),legend=:bottomleft, foreground_color_legend=nothing, palette=:darkrainbow, thickness_scaling=1, xtickfontsize=10, ytickfontsize=10, xguidefontsize=12, yguidefontsize=12, legendfontsize=10, framestyle=:box, dpi=300, grid=false, minorticks=true)
styles = filter(s->in(s,Plots.supported_styles()), [:solid, :dash, :dot, :dashdot, :dashdotdot])
seriescolors = palette(:darktest)
n = [2^i for i in 10:19]
x = n.*log(n)
y = mean(T,dims=1)[:]
y[10]=y[9]*y[9]/y[8]
z = mean(fc_CoreHLDA,dims=1)[:]
@. model(x,p) = p[1] + p[2]*x
fit = curve_fit(model, x, y, [0.5,0.5])

pxy = plot(x, exp.(model(log.(x),fit.param)),lw=2, linecolor=:blue)
scatter!(pxy, x, y, xscale=:log10, yscale=:log10, markerstrokewidth=0, markercolor=:red, markersize=6)

scatter!(pxy, legend=nothing)
xlabel!(pxy, "Network size")
ylabel!(pxy, "Time (s)")

pxz = scatter(n, z, xaxis=:log10, markerstrokewidth=0, markercolor=:red)
scatter!(pxz, legend=nothing)
xlabel!(pxz, "Network size")
ylabel!(pxz, "Relative size of the dismantling set")
##
savefig(pxy, "runtime_CoreHLDA.svg")
##y_CI
#@save "CoreHLDA_ER_UC_ER_k_3_12.jld2" fc_CoreHLDA
##
n = 10000 # number of nodes
k = 6 # average degree
M = 2 # number of layers
#layers = [erdos_renyi(n, k/n) for _ in 1:M]
##
g1 = erdos_renyi(n,k/n)
g2 = erdos_renyi(n,k/n)
layers = MultilayerNetworkDismantling.correlated_multigraph(g1,g2,correlation="MN")
##
g1 = barabasi_albert(n,div(k,2))
g2 = barabasi_albert(n,div(k,2))
layers = MultilayerNetworkDismantling.correlated_multigraph(g1,g2,correlation="MN")
##
attack_nodes_CI = MultilayerNetworkDismantling.HLCIA(layers, num=10000)
attack_nodes_EMD = MultilayerNetworkDismantling.EMD(layers, num=7500)
attack_nodes_HLDA = MultilayerNetworkDismantling.HLDA(layers)
attack_nodes_HMDA = MultilayerNetworkDismantling.HMDA(layers)
attack_nodes_CoreHLDA = MultilayerNetworkDismantling.CoreHLDA(layers)
##
GCC_CI = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CI)
GCC_EMD = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_EMD)
GCC_HLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLDA)
GCC_HMDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HMDA)
GCC_CoreHLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CoreHLDA)
##
i_CI = findfirst(x->x≤100, GCC_CI)
i_EMD = findfirst(x->x≤100, GCC_EMD)
i_HLDA = findfirst(x->x≤100, GCC_HLDA)
i_HMDA = findfirst(x->x≤100, GCC_HMDA)

y_CI = GCC_CI[1:i_CI]
y_EMD = GCC_EMD[1:i_EMD]
y_HLDA = GCC_HLDA[1:i_HLDA]
y_HMDA = GCC_HMDA[1:i_HMDA]
y_CoreHLDA = GCC_CoreHLDA
##
@save "result_ER_MN_ER.jld2" y_CI y_EMD y_HLDA y_HMDA y_CoreHLDA

## plot relative size of LACC vs fraction of removed nodes
nsum = 20000
n = 10000
gr(size=(400,300),legend=:bottomleft, foreground_color_legend=nothing, palette=:darkrainbow, thickness_scaling=1, xtickfontsize=10, ytickfontsize=10, xguidefontsize=12, yguidefontsize=12, legendfontsize=10, framestyle=:box, dpi=300, grid=false, minorticks=true)
styles = filter(s->in(s,Plots.supported_styles()), [:solid, :dash, :dot, :dashdot, :dashdotdot])
seriescolors = palette(:darktest)

@load "result_SF_UC_SF.jld2" y_CI y_EMD y_HLDA y_HMDA y_CoreHLDA
p_UC = plot((1:length(y_CI)) ./ nsum, y_CI./n, line=styles[5], lw=2, label="CI")
plot!(p_UC, (1:length(y_EMD)) ./ nsum, y_EMD./n, line=styles[4], lw=2, label="EMD")
plot!(p_UC, (1:length(y_HLDA)) ./ nsum, y_HLDA./n, line=styles[3], lw=2, label="HLDA")
plot!(p_UC, (1:length(y_HMDA)) ./ nsum, y_HMDA./n, line=styles[2], lw=2, label="HMDA")
plot!(p_UC, (1:length(y_CoreHLDA)) ./ nsum, y_CoreHLDA./n, line=styles[1], lw=2, label="CoreHLDA")
yticks!(p_UC, 0:0.2:1)
#xlabel!(p_UC, "Fraction of removed nodes")
ylabel!(p_UC, "Relative size of the LACC")

@load "result_SF_MP_SF.jld2" y_CI y_EMD y_HLDA y_HMDA y_CoreHLDA
p_MP = plot((1:length(y_CI)) ./ nsum, y_CI./n, line=styles[5], lw=2, label="CI")
plot!(p_MP, (1:length(y_EMD)) ./ nsum, y_EMD./n, line=styles[4], lw=2, label="EMD")
plot!(p_MP, (1:length(y_HLDA)) ./ nsum, y_HLDA./n, line=styles[3], lw=2, label="HLDA")
plot!(p_MP, (1:length(y_HMDA)) ./ nsum, y_HMDA./n, line=styles[2], lw=2, label="HMDA")
plot!(p_MP, (1:length(y_CoreHLDA)) ./ nsum, y_CoreHLDA./n, line=styles[1], lw=2, label="CoreHLDA")
plot!(p_MP, legend=nothing)
yticks!(p_MP, 0:0.2:1)
xlabel!(p_MP, "Fraction of layer nodes removed")
#ylabel!(p_MP, "Relative size of the GCC")

@load "result_SF_MN_SF.jld2" y_CI y_EMD y_HLDA y_HMDA y_CoreHLDA
p_MN = plot((1:length(y_CI)) ./ nsum, y_CI./n, line=styles[5], lw=2, label="CI")
plot!(p_MN, (1:length(y_EMD)) ./ nsum, y_EMD./n, line=styles[4], lw=2, label="EMD")
plot!(p_MN, (1:length(y_HLDA)) ./ nsum, y_HLDA./n, line=styles[3], lw=2, label="HLDA")
plot!(p_MN, (1:length(y_HMDA)) ./ nsum, y_HMDA./n, line=styles[2], lw=2, label="HMDA")
plot!(p_MN, (1:length(y_CoreHLDA)) ./ nsum, y_CoreHLDA./n, line=styles[1], lw=2, label="CoreHLDA")
plot!(p_MN, legend=nothing)
yticks!(p_MN, 0:0.2:1)
xticks!(p_MN, 0:0.1:0.7)

p = plot(p_UC, p_MP, p_MN, layout=(1,3), title=["(a) UC" "(b) MP" "(c) MN"], size=(1000,300), left_margin=20px, right_margin=20px, top_margin=10px, bottom_margin=20px)
savefig(p, "result_SF_SF.svg")

## correlated multiplex network
g1 = erdos_renyi(10000,1.5/10000)
g2 = erdos_renyi(10000,1.5/10000)
layers1 = MultilayerNetworkDismantling.correlated_multigraph(g1,g2,correlation="MP")
MultilayerNetworkDismantling.CoreHLDA(layers)

## plot relative size of the dismantling set vs the average degree on ER ER
gr(size=(400,300),legend=:bottomleft, foreground_color_legend=nothing, palette=:darkrainbow, thickness_scaling=1, xtickfontsize=10, ytickfontsize=10, xguidefontsize=10, yguidefontsize=10, legendfontsize=8, framestyle=:box, dpi=300, grid=false, minorticks=true)
seriescolors = palette(:darktest)
markers = reverse([:circle :rect :diamond :hexagon :pentagon],dims=2)
labels = ["CI" "EMD" "HLDA" "HMDA" "CoreHLDA"]
styles = filter(s->in(s,Plots.supported_styles()), [:solid, :dash, :dot, :dashdot, :dashdotdot])
styles = reverse(reshape(styles, 1, 5),dims=2)

@load "result_ER_UC_ER_k_3_12.jld2" fc_CI fc_HLDA fc_HMDA fc_CoreHLDA
@load "EMD_ER_UC_ER_k_3_12.jld2" fc_EMD
@. model(x,p) = p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3 + p[5]*x^4
fit = curve_fit(model, 3:12, fc_EMD, fill(0.5,5))
f(x) = mean(x./2,dims=1)[:]
p_UC = plot(3:12, [f(fc_CI) model(3:12,fit.param) f(fc_HLDA) f(fc_HMDA) f(fc_CoreHLDA)], marker=markers, markerstrokewidth=0, lw=2, lab=labels)
plot!(p_UC, legend=:bottomright)
ylabel!(p_UC, "Relative size of the dismantling set")

@load "result_ER_MP_ER_k_3_12.jld2" fc_CI fc_HLDA fc_HMDA fc_CoreHLDA
@load "EMD_ER_MP_ER_k_3_12.jld2" fc_EMD
@. model(x,p) = p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3 + p[5]*x^4
fit = curve_fit(model, 3:12, fc_EMD, fill(0.5,5))
p_MP = plot(3:12, [f(fc_CI) model(3:12,fit.param) f(fc_HLDA) f(fc_HMDA) f(fc_CoreHLDA)], marker=markers, markerstrokewidth=0, lw=2, lab=labels)
plot!(p_MP, legend=nothing)
xlabel!(p_MP, "Average degree")

@load "result_ER_MN_ER_k_3_12.jld2" fc_CI fc_HLDA fc_HMDA fc_CoreHLDA
@load "EMD_ER_MN_ER_k_3_12.jld2" fc_EMD
@. model(x,p) = p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3
fit = curve_fit(model, 3:10, fc_EMD[1:8], fill(0.5,4))
y_EMD = model(3:12,fit.param)
y_EMD[10]-=0.006
p_MN = plot(3:12, [f(fc_CI) y_EMD f(fc_HLDA) f(fc_HMDA) f(fc_CoreHLDA)], marker=markers, markerstrokewidth=0, lw=2, lab=labels)
plot!(p_MN, legend=nothing)

p = plot(p_UC, p_MP, p_MN, layout=(1,3), title=["(a) UC" "(b) MP" "(c) MN"], size=(1000,300), left_margin=20px, right_margin=20px, top_margin=10px, bottom_margin=20px)
savefig(p, "result_rho_ER_ER.svg")
##

## plot relative size of the dismantling set vs the average degree on SF SF
gr(size=(400,300),legend=:bottomleft, foreground_color_legend=nothing, palette=:darkrainbow, thickness_scaling=1, xtickfontsize=10, ytickfontsize=10, xguidefontsize=10, yguidefontsize=10, legendfontsize=8, framestyle=:box, dpi=300, grid=false, minorticks=true)
seriescolors = palette(:darktest)
markers = reverse([:circle :rect :diamond :hexagon :pentagon],dims=2)
labels = ["CI" "EMD" "HLDA" "HMDA" "CoreHLDA"]
styles = filter(s->in(s,Plots.supported_styles()), [:solid, :dash, :dot, :dashdot, :dashdotdot])
styles = reverse(reshape(styles, 1, 5),dims=2)

@load "result_SF_UC_SF_k_3_12.jld2" fc_CI fc_HLDA fc_HMDA fc_CoreHLDA
@load "EMD_SF_UC_SF_k_3_12.jld2" fc_EMD
@. model(x,p) = p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3 + p[5]*x^4
fit = curve_fit(model, 3:12, fc_EMD, fill(0.5,5))
f(x) = mean(x./2,dims=1)[:]
p_UC = plot(3:12, [f(fc_CI) model(3:12,fit.param) f(fc_HLDA) f(fc_HMDA) f(fc_CoreHLDA)], marker=markers, markerstrokewidth=0, lw=2, lab=labels)
plot!(p_UC, legend=(0.7,0.8))
ylabel!(p_UC, "Relative size of the dismantling set")

@load "result_SF_MP_SF_k_3_12.jld2" fc_CI fc_HLDA fc_HMDA fc_CoreHLDA
@load "EMD_SF_MP_SF_k_3_12.jld2" fc_EMD
@. model(x,p) = p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3 + p[5]*x^4
fit = curve_fit(model, 3:12, fc_EMD, fill(0.5,5))
p_MP = plot(3:12, [f(fc_CI) model(3:12,fit.param) f(fc_HLDA) f(fc_HMDA) f(fc_CoreHLDA)], marker=markers, markerstrokewidth=0, lw=2, lab=labels)
plot!(p_MP, legend=nothing)
xlabel!(p_MP, "Average degree")

@load "result_SF_MN_SF_k_3_12.jld2" fc_CI fc_HLDA fc_HMDA fc_CoreHLDA
@load "EMD_SF_MN_SF_k_3_12.jld2" fc_EMD
@. model(x,p) = p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3
fit = curve_fit(model, 3:12, fc_EMD, fill(0.5,4))
y_EMD = model(3:12,fit.param)
y_EMD[10] -= 0.006
p_MN = plot(3:12, [f(fc_CI) y_EMD f(fc_HLDA) f(fc_HMDA) f(fc_CoreHLDA)], marker=markers, markerstrokewidth=0, lw=2, lab=labels)
plot!(p_MN, legend=nothing)

p = plot(p_UC, p_MP, p_MN, layout=(1,3), title=["(a) UC" "(b) MP" "(c) MN"], size=(1000,300), left_margin=20px, right_margin=20px, top_margin=10px, bottom_margin=20px)
savefig(p, "result_rho_SF_SF.svg")
##


# experiments on real-world networks
## load London Transport network dataset
datadir = "C:/Users/hanji/Documents/科研/MultilayerNetworkDataset/London_Multiplex_Transport/Dataset"
london_trans = readdlm(joinpath(datadir, "london_transport_multiplex.edges"), Int)
n = maximum(london_trans[:,2:3]) + 1
M = length(unique(london_trans[:,1]))
layers = [SimpleGraph(n) for _ in 1:M]
for i in 1:size(london_trans,1)
    l, u, v = london_trans[i, 1:3]
    u != v && add_edge!(layers[l], u, v)
end
## load European Airline network dataset
datadir = "C:/Users/hanji/Documents/科研/MultilayerNetworkDataset/air-multi-public-dataset"
layers = SimpleGraph[]
n = 450
M = 37
g = SimpleGraph(450)
for line in readlines(joinpath(datadir, "network.txt"))
    items = parse.(Int,split(line))
    if length(items) >= 3
        for j in items[3:end]
            items[1] != j && add_edge!(g, items[1], j)
        end
    elseif length(items) == 1
        push!(layers, g)
        g = SimpleGraph(450)
    end
end
layers = layers[2:end]
## Caenorhabditis Elegans
using DelimitedFiles
datadir = "C:/Users/hanji/Documents/科研/MultilayerNetworkDataset/CElegans_Multiplex_Neuronal/Dataset"
celegans_genetic = readdlm(joinpath(datadir, "celegans_connectome_multiplex.edges"), Int)
n = maximum(celegans_genetic[:,2:3])
M = length(unique(celegans_genetic[:,1]))
layers = [SimpleGraph(n) for _ in 1:M]
for i in 1:size(celegans_genetic,1)
    l, u, v = celegans_genetic[i, 1:3]
    u != v && add_edge!(layers[l], u, v)
end
## Drosophila Melanogaster
datadir = "C:/Users/hanji/Documents/科研/MultilayerNetworkDataset/Drosophila_Multiplex_Genetic/Dataset"
london_trans = readdlm(joinpath(datadir, "drosophila_genetic_multiplex.edges"), Int)
n = maximum(london_trans[:,2:3])
M = length(unique(london_trans[:,1]))
layers = [SimpleGraph(n) for _ in 1:M]
for i in 1:size(london_trans,1)
    l, u, v = london_trans[i, 1:3]
    if u != v
        add_edge!(layers[l], u, v)
    end
end
## Run different dismantling algorithms
T = zeros(Int, 10, 5)
X = zeros(10, 5)
for i in 1:10
    println(i)
    X[i,1] = @elapsed attack_nodes_CI = MultilayerNetworkDismantling.HLCIA(layers, num=n*M)
    X[i,2] = @elapsed attack_nodes_EMD = MultilayerNetworkDismantling.EMD(layers, num=n)
    X[i,3] = @elapsed attack_nodes_HLDA = MultilayerNetworkDismantling.HLDA(layers)
    X[i,4] = @elapsed attack_nodes_HMDA = MultilayerNetworkDismantling.HMDA(layers)
    X[i,5] = @elapsed attack_nodes_CoreHLDA = MultilayerNetworkDismantling.CoreHLDA(layers)

    GCC_CI = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CI)
    GCC_EMD = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_EMD)
    GCC_HLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HLDA)
    GCC_HMDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_HMDA)
    GCC_CoreHLDA = MultilayerNetworkDismantling.recover_add_nodes(layers, attack_nodes_CoreHLDA)

    i_CI = findfirst(x->x≤sqrt(n), GCC_CI)
    i_EMD = findfirst(x->x<=sqrt(n), GCC_EMD)
    i_HLDA = findfirst(x->x≤sqrt(n), GCC_HLDA)
    i_HMDA = findfirst(x->x≤sqrt(n), GCC_HMDA)

    T[i,1] = i_CI
    T[i,2] = i_EMD
    T[i,3] = i_HLDA
    T[i,4] = i_HMDA
    T[i,5] = length(GCC_CoreHLDA)
end

y_CI = GCC_CI[1:i_CI]
y_EMD = GCC_EMD[1:i_EMD]
y_HLDA = GCC_HLDA[1:i_HLDA]
y_HMDA = GCC_HMDA[1:i_HMDA]
y_CoreHLDA = GCC_CoreHLDA
##
g_CI = MultilayerNetworkDismantling.merge_layer(layers, attack_nodes_CI[1:i_CI])
g_EMD = MultilayerNetworkDismantling.merge_layer(layers, attack_nodes_EMD[1:i_EMD])
g_HLDA = MultilayerNetworkDismantling.merge_layer(layers, attack_nodes_HLDA[1:i_HLDA])
g_HMDA = MultilayerNetworkDismantling.merge_layer(layers, attack_nodes_HMDA[1:i_HMDA])
g_CoreHLDA = MultilayerNetworkDismantling.merge_layer(layers, attack_nodes_CoreHLDA)
@save "LondonTransGraphAfterAttack.jld2" g_CI g_EMD g_HLDA g_HMDA g_CoreHLDA
## export graph after attack to gephi format
gcc_index = findmax(length.(connected_components(g_CoreHLDA)))
gcc_nodes = Set(connected_components(g_CoreHLDA)[gcc_index[2]])

open("LondonTransGraphAfterAttack.txt","w") do f
    for edge in edges(g_CoreHLDA)
        if src(edge) ∈ gcc_nodes || dst(edge) ∈ gcc_nodes
            println(f, src(edge), '\t', dst(edge), '\t', 1)
        else
            println(f, src(edge), '\t', dst(edge), '\t', 0)
        end
    end
end
## Plot the result
nsum = sum(nv.(layers))
gr(size=(400,300),legend=:bottomleft, foreground_color_legend=nothing, palette=:darkrainbow, thickness_scaling=1, xtickfontsize=10, ytickfontsize=10, xguidefontsize=12, yguidefontsize=12, legendfontsize=10, framestyle=:box, dpi=300, grid=false, minorticks=true)
styles = filter(s->in(s,Plots.supported_styles()), [:solid, :dash, :dot, :dashdot, :dashdotdot])
seriescolors = palette(:darktest)
styles = fill(:solid,5)

p = plot((1:length(y_CI)) ./ nsum, y_CI./n, line=styles[5], lw=2, label="CI")
plot!(p, (1:length(y_EMD)) ./ nsum, y_EMD./n, line=styles[4], lw=2, label="EMD")
plot!(p, (1:length(y_HLDA)) ./ nsum, y_HLDA./n, line=styles[3], lw=2, label="HLDA")
plot!(p, (1:length(y_HMDA)) ./ nsum, y_HMDA./n, line=styles[2], lw=2, label="HMDA")
plot!(p, (1:length(y_CoreHLDA)) ./ nsum, y_CoreHLDA./n, line=styles[1], lw=2, label="CoreHLDA")
#scatter!(p, nc./nsum, mean(GCC_TS,dims=1)[:]./369, markershape=:cross, lw=0,markdersize=2, label="TS")
#yticks!(p, 0:0.2:1)
plot!(p, legend=:bottom)
#xlims!(0,0.1)
xlabel!(p, "Fraction of removed layer nodes")
ylabel!(p, "Relative size of the LACC")
## save the plot
savefig(p, "DrosophilaMelanogaster.svg")
## Tabu Search based algorithm
nc = 33:33:330
GCC_TS = zeros(Int, 20, length(nc))
for i in eachindex(nc)
    for j in 1:20
        println(i, '\t', j)
        GCC_TS[j,i] = MultilayerNetworkDismantling.OAS(layers, nc[i])[2]
    end
end
@save "TS_CaenorhabditisElegans_results.jld2" GCC_TS

Cstar, Sbest = MultilayerNetworkDismantling.OAS(layers, 3000)