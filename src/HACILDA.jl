function union_degree(layers, presents, i, mask, neis)
    nei = Set{Int}()
    index = 0
	for l in eachindex(layers)
		if presents[l][i]
			for j in neighbors(layers[l],i)
                if presents[l][j] && !mask[j]
                    mask[j] = true
                    index += 1
                    neis[index] = j
				end
			end
		end
    end
    for k in 1:index
        mask[neis[k]] = false
    end
	return index
end

"""
compute collective influence score
"""
function collective_influence1!(g, u, scope, bfsqueue, is_visited, present, deg)
	deg[u] == 0 && return 0
	scope == 0 && return deg[u]
	_, startIt, endIt = frontier_scope1!(g, u, scope, bfsqueue, is_visited, present)
	ci = 0
	for i=startIt:endIt
		ci += deg[bfsqueue[i]] - 1
	end
	ci *= deg[u] - 1
	return ci
end

"""
fontier is between startIt and endIt
"""
function frontier_scope1!(g, u, scope, bfsqueue, is_visited, present)
    startIt1, startIt, endIt  = 1, 1, 2
	bfsqueue[1] = u
	is_visited[u] = true

	for i=1:scope
		lastEndIt = endIt
		while startIt != lastEndIt
			for v in neighbors(g, bfsqueue[startIt])
				if present[v] && !is_visited[v]
					bfsqueue[endIt] = v
					endIt += 1
					is_visited[v] = true
				end
			end
			startIt += 1
        end
        if i == 1
            startIt1 = endIt-1
        end
	end

	for i=1:endIt-1
		is_visited[bfsqueue[i]] = false
	end
	return startIt1, startIt, endIt-1
end

function HACILDA(layers, T=sum(nv.(layers)), scope=1)
    g = foldl(union, layers)
    n = nv(g)
    mask = fill(false, n)
    neis = zeros(Int, n)
	deg = degree(g)
    attacked_nodes = @NamedTuple{layer::Int,id::Int}[]
	presents = [fill(true, n) for _ in eachindex(layers)]
    degs = degree.(layers)

    bfsqueue = zeros(Int, n)
	bfsqueue1 = zeros(Int, n)
    is_visited = fill(false, n)
    present = fill(true, n)
    ci = [collective_influence1!(g, u, scope, bfsqueue, is_visited, present, deg) for u in vertices(g)]

	H = SortedDict{Tuple{Int,Int},Set{eltype(attacked_nodes)}}()

	for l in eachindex(layers)
		for u in vertices(layers[l])
			score_u = (ci[u],degs[l][u])
			if haskey(H, score_u)
				push!(H[score_u], (layer=l,id=u))
			else
				push!(H, score_u=>Set([(layer=l,id=u)]))
			end
		end
	end

	score_max = last(H)[1]
	done = false
	for t in 1:T
		i = rand(H[score_max])
		delete!(H[score_max],i)
        presents[i.layer][i.id] = false
        present[i.id] = false
        for l in eachindex(layers)
            if presents[l][i.id]
                present[i.id] = true
                break
            end
        end
        push!(attacked_nodes, i)
        
        new_deg_i = union_degree(layers, presents, i.id, mask, neis)
        if new_deg_i < deg[i.id]
            startIt1, startIt, endIt = frontier_scope1!(g, i.id, scope+1, bfsqueue1, is_visited, present)
            for index in 1:endIt
                j = bfsqueue1[index]
                for l in eachindex(layers)
                    if presents[l][j]
                        score_j = (ci[j], degs[l][j])
                        delete!(H[score_j], (layer=l, id=j))
                    end
                end
            end

            deg[i.id] = new_deg_i
            for j in neighbors(layers[i.layer], i.id)
                if presents[i.layer][j]
                    degs[i.layer][j] -= 1
                    deg[j] = union_degree(layers, presents, j, mask, neis)
                end
            end

            for index in 1:endIt
                j = bfsqueue1[index]
                ci[j] = collective_influence1!(g, j, scope, bfsqueue, is_visited, present, deg)
                for l in eachindex(layers)
                    if presents[l][j]
                        score_j = (ci[j], degs[l][j])
                        if haskey(H, score_j)
                            push!(H[score_j], (layer=l,id=j))
                        else
                            push!(H, score_j=>Set([(layer=l,id=j)]))
                        end
                    end
                end
            end
        else
            for j in neighbors(layers[i.layer], i.id)
                if presents[i.layer][j]
                    score_j = (ci[j], degs[i.layer][j])
                    delete!(H[score_j], (layer=i.layer, id=j))
                    degs[i.layer][j] -= 1
                    score_j = (ci[j], degs[i.layer][j])
                    if haskey(H, score_j)
                        push!(H[score_j], (layer=i.layer,id=j))
                    else
                        push!(H, score_j=>Set([(layer=i.layer,id=j)]))
                    end
                end
            end
        end
        #println(score_max)
		while isempty(H[score_max])
            delete!(H,score_max)
            if isempty(H)
				done = true
				break
			end
            score_max = last(H)[1]
            if score_max ≤ (0,0)
                done = true
                break
            end
		end
		done && break
	end
	attacked_nodes
end