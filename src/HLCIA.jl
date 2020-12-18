"""
compute collective influence score
"""
function collective_influence!(g, u, scope, bfsqueue, is_visited, present, deg)
	deg[u] == 0 && return 0
	scope == 0 && return deg[u]
	startIt, endIt = frontier_scope!(g, u, scope, bfsqueue, is_visited, present)
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
function frontier_scope!(g, u, scope, bfsqueue, is_visited, present)
    startIt, endIt = 1, 2
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
	end

	for i=1:endIt-1
		is_visited[bfsqueue[i]] = false
	end
	return startIt, endIt-1
end

function HLCIA(layers; scope=1, num=round(Int,sqrt(nv(layers[1]))))
    g = foldl(union, layers)
	n = nv(g)
    attack_nodes = @NamedTuple{layer::Int,id::Int}[]
	bfsqueue = zeros(Int, n)
	bfsqueue1 = zeros(Int, n)
	is_visited = fill(false, n)
	present = [fill(true, n) for _ in eachindex(layers)]
    deg = degree.(layers)
	ci = [zeros(Int, n) for _ in eachindex(layers)]
    # a sorted dict with CI score as keys and
    # sets of nodes with the same CI as values
    sdict = SortedDict{Int,Set{eltype(attack_nodes)}}()
    for l in eachindex(layers)
        for u in vertices(layers[l])
            ci[l][u] = collective_influence!(layers[l], u, scope, bfsqueue, is_visited, present[l], deg[l])
            if haskey(sdict, ci[l][u])
                push!(sdict[ci[l][u]], (layer=l,id=u))
            else
                push!(sdict, ci[l][u]=>Set([(layer=l,id=u)]))
            end
        end
    end
	t = 1
	while !isempty(sdict) && t <= num
		k, candidates = last(sdict)
		if !isempty(candidates)
			t += 1
            u = rand(candidates)
            push!(attack_nodes, u)
            
			startIt, endIt = frontier_scope!(layers[u.layer], u.id, scope+1, bfsqueue, is_visited, present[u.layer])
			delete!(candidates,u)

			# remove node u
			present[u.layer][u.id] = false
			for v in neighbors(layers[u.layer],u.id)
				if present[u.layer][v]
					deg[u.layer][v] -= 1
				end
			end

			# update CI score of neighbors and update sdict
			for i in startIt:endIt
				v = bfsqueue[i]
				delete!(sdict[ci[u.layer][v]], (layer=u.layer,id=v))
				ci[u.layer][v] = ci[u.layer][v] - deg[u.layer][v] + 1
				if haskey(sdict, ci[u.layer][v])
					push!(sdict[ci[u.layer][v]], (layer=u.layer,id=v))
				else
					push!(sdict, ci[u.layer][v]=>Set([(layer=u.layer,id=v)]))
				end
			end
			for i in 2:startIt-1
				v = bfsqueue[i]
				delete!(sdict[ci[u.layer][v]], (layer=u.layer,id=v))
				ci[u.layer][v] = collective_influence!(layers[u.layer], v, scope, bfsqueue1, is_visited, present[u.layer], deg[u.layer])
				if haskey(sdict, ci[u.layer][v])
					push!(sdict[ci[u.layer][v]], (layer=u.layer,id=v))
				else
					push!(sdict, ci[u.layer][v]=>Set([(layer=u.layer,id=v)]))
				end
			end
		else
			pop!(sdict, k)
		end
    end
    attack_nodes
end