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

function HMCIA(layers; scope=1, num=round(Int,sqrt(nv(layers[1]))))
    g = foldl(union, layers)
	n = nv(g)
    attack_nodes = @NamedTuple{layer::Int,id::Int}[]
	bfsqueue = zeros(Int, n)
	bfsqueue1 = zeros(Int, n)
	is_visited = fill(false, n)
	present = fill(true, n)
    deg = degree(g)
	ci = zeros(Int, n)
    # a sorted dict with CI score as keys and
    # sets of nodes with the same CI as values
    sdict = SortedDict{Int,Set{Int}}()
    for u in vertices(g)
		ci[u] = collective_influence!(g, u, scope, bfsqueue, is_visited, present, deg)
		if haskey(sdict, ci[u])
            push!(sdict[ci[u]], u)
        else
            push!(sdict, ci[u]=>Set(u))
        end
    end
	t = 1
	while !isempty(sdict) && t <= num
		k, candidates = last(sdict)
		if !isempty(candidates)
			t += 1
            u = rand(candidates)
            for l in eachindex(layers)
                push!(attack_nodes, (layer=l,id=u))
            end
			startIt, endIt = frontier_scope!(g, u, scope+1, bfsqueue, is_visited, present)
			delete!(candidates,u)

			# remove node u
			present[u] = false
			for v in neighbors(g,u)
				if present[v]
					deg[v] -= 1
				end
			end

			# update CI score of neighbors and update sdict
			for i in startIt:endIt
				v = bfsqueue[i]
				delete!(sdict[ci[v]], v)
				ci[v] = ci[v] - deg[v] + 1
				if haskey(sdict, ci[v])
					push!(sdict[ci[v]], v)
				else
					push!(sdict, ci[v]=>Set(v))
				end
			end
			for i in 2:startIt-1
				v = bfsqueue[i]
				delete!(sdict[ci[v]], v)
				ci[v] = collective_influence!(g, v, scope, bfsqueue1, is_visited, present, deg)
				if haskey(sdict, ci[v])
					push!(sdict[ci[v]], v)
				else
					push!(sdict, ci[v]=>Set(v))
				end
			end
		else
			pop!(sdict, k)
		end
    end
    attack_nodes
end