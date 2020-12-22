function OAS(layers, n, L=10, nc=100, maxiter=1000)
    S = zeros(Int, nc)
    ix = Vector{Int}(undef, nc)
    label = zeros(Int,nv(layers[1]))
	cnt = zeros(Int,nv(layers[1]))
    nodes = [(l,i) for l in eachindex(layers), i in vertices(layers[1])]
    C0 = [fill(true,nv(layers[1])) for _ in eachindex(layers)]
    attack_nodes = Set(sample(nodes,n,replace=false))
    present_nodes = Set(setdiff(nodes, attack_nodes))
    for (l,i) in attack_nodes
        C0[l][i] = false
    end
    tabulist = Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}}, Int}()
    Cstar = deepcopy(C0)
    Cnow = deepcopy(C0)
    Sbest = Snow = MLCC!(layers, C0, label, cnt)
    iter = 0
    while iter < maxiter
        neilist = Tuple{Tuple{Int,Int},Tuple{Int,Int}}[]
        for num in 1:nc
            node1 = rand(attack_nodes)
            node2 = rand(present_nodes)
            push!(neilist, (node1,node2))
            Cnow[node1[1]][node1[2]], Cnow[node2[1]][node2[2]] = Cnow[node2[1]][node2[2]], Cnow[node1[1]][node1[2]]
            S[num] = MLCC!(layers, Cnow, label, cnt)
            Cnow[node1[1]][node1[2]], Cnow[node2[1]][node2[2]] = Cnow[node2[1]][node2[2]], Cnow[node1[1]][node1[2]]
        end
        index = partialsortperm!(ix, S, 1:2)
        node1, node2 = neilist[index[1]]
        if !haskey(tabulist, neilist[index[1]]) || S[index[1]] < Snow
            node1, node2 = neilist[index[1]]
            Cnow[node1[1]][node1[2]], Cnow[node2[1]][node2[2]] = Cnow[node2[1]][node2[2]], Cnow[node1[1]][node1[2]]
            Snow = S[index[1]]
            push!(tabulist, (node1,node2)=>L)
            delete!(attack_nodes, node1)
            push!(present_nodes, node1)
            delete!(present_nodes, node2)
            push!(attack_nodes, node2)
        elseif !haskey(tabulist, neilist[index[2]]) || S[index[2]] < Snow
            node1, node2 = neilist[index[2]]
            Cnow[node1[1]][node1[2]], Cnow[node2[1]][node2[2]] = Cnow[node2[1]][node2[2]], Cnow[node1[1]][node1[2]]
            Snow = S[index[2]]
            push!(tabulist, (node1,node2)=>L)
            delete!(attack_nodes, node1)
            push!(present_nodes, node1)
            delete!(present_nodes, node2)
            push!(attack_nodes, node2)
        else
            continue
        end
        if Snow < Sbest
            Cstar[node1[1]][node1[2]], Cstar[node2[1]][node2[2]] = Cstar[node2[1]][node2[2]], Cstar[node1[1]][node1[2]]
            Sbest = Snow
            iter = 0
        else
            iter += 1
        end
        for k in keys(tabulist)
            tabulist[k] -= 1
            if tabulist[k] < 1
                delete!(tabulist, k)
            end
        end
    end
    Cstar, Sbest
end