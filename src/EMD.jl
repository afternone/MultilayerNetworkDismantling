function effective_multiplex_degree!(w, tempw, g, layers, degs, deg, present, T=10)
    M = length(layers)
    for i in vertices(g)
        if present[i]
            tempw[i] = w[i] = deg[i]
        else
            tempw[i] = w[i] = 0
        end
    end
    for t in 1:T
        for i in vertices(g)
            if present[i]
                wi = 0
                for l in eachindex(layers)
                    for j in neighbors(layers[l],i)
                        if present[j] && degs[l][j] > 0
                            wi += w[j] / degs[l][j]
                        end
                    end
                end
                tempw[i] = wi / M
            end
        end
        w[:] = tempw[:]
    end
end

function EMD(layers; num=nv(layers[1]), T=10)
    attack_nodes = @NamedTuple{layer::Int,id::Int}[]
    g = foldl(union, layers)
    deg = degree(g)
    degs = degree.(layers)
    w = zeros(nv(g))
    tempw = zeros(nv(g))
    present = fill(true, nv(g))

    for _ in 1:num
        effective_multiplex_degree!(w, tempw, g, layers, degs, deg, present, T)
        wmax = -1.0
        ibest = 0
        for i in vertices(g)
            if present[i] && w[i] > wmax
                wmax = w[i]
                ibest = i
            end
        end
        for l in eachindex(layers)
            push!(attack_nodes, (layer=l, id=ibest))
            for j in neighbors(layers[l], ibest)
                if present[j]
                    degs[l][j] -= 1
                end
            end
        end
        for j in neighbors(g,ibest)
            if present[j]
                deg[j] -= 1
            end
        end
        present[ibest] = false
    end
	return attack_nodes
end
