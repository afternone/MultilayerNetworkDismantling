function MLCC!(layers, presents, label = zeros(Int,nv(layers[1])), cnt = zeros(Int,nv(layers[1])))
    n = nv(layers[1])
    Q = Queue{Int}()
	maxcnt = 0
	for u in vertices(layers[1])
		is_present = false
		for l in eachindex(layers)
			if presents[l][u]
				is_present = true
				break
			end
		end
		if is_present
			label[u] != 0 && continue
			label[u] = u
			cnt[u] += 1
			enqueue!(Q,u)
			while !isempty(Q)
				src_node = dequeue!(Q)
				for l in eachindex(layers)
					if presents[l][src_node]
						for v in all_neighbors(layers[l],src_node)
							if presents[l][v] && label[v] == 0
								enqueue!(Q, v)
								label[v] = u
								cnt[u] += 1
							end
						end
					end
				end
			end
			if cnt[u] > maxcnt
				maxcnt = cnt[u]
			end
		end
	end
	for i in eachindex(label)
		label[i] = 0
		cnt[i] = 0
	end
	return maxcnt
end

function LCC(g, present)
    n = nv(g)
    Q = Queue{Int}()
	label = zeros(Int,n)
	cnt = zeros(Int,n)
	maxcnt = 0
    @inbounds for u in vertices(g)
		if present[u]
	        label[u] != 0 && continue
	        label[u] = u
			cnt[u] += 1
	        enqueue!(Q, u)
	        while !isempty(Q)
	            src = dequeue!(Q)
	            for vertex in all_neighbors(g, src)
					if present[vertex]
		                if label[vertex] == 0
		                    enqueue!(Q, vertex)
		                    label[vertex] = u
							cnt[u] += 1
		                end
					end
	            end
	        end
			if cnt[u] > maxcnt
				maxcnt = cnt[u]
			end
		end
    end
    return maxcnt
end

function LCC(g)
    n = nv(g)
    Q = Queue{Int}()
	label = zeros(Int,n)
	cnt = zeros(Int,n)
	maxcnt = 0
    @inbounds for u in vertices(g)
		label[u] != 0 && continue
		label[u] = u
		cnt[u] += 1
		enqueue!(Q, u)
		while !isempty(Q)
			src = dequeue!(Q)
			for vertex in all_neighbors(g, src)
				if label[vertex] == 0
					enqueue!(Q, vertex)
					label[vertex] = u
					cnt[u] += 1
				end
			end
		end
		if cnt[u] > maxcnt
			maxcnt = cnt[u]
		end
    end
    return maxcnt
end

function merge_layer(layers, presents::Array{Array{Bool,1},1})
	g = SimpleGraph(nv(layers[1]))
	for (layer, present) in zip(layers, presents)
		for e in edges(layer)
			if !has_edge(g, e) && present[src(e)] && present[dst(e)]
				add_edge!(g, e)
			end
		end
	end
	return g
end

function merge_layer(layers, attack_nodes::Array{NamedTuple{(:layer, :id),Tuple{Int64,Int64}},1})
	presents = [fill(true, nv(layers[1])) for _ in eachindex(layers)]
	for i in attack_nodes
		presents[i.layer][i.id] = false
	end
	merge_layer(layers, presents)
end

function merge_nodes!(ds, comp_sizes, i, j)
	iroot = find_root(ds, i)
	jroot = find_root(ds, j)
	if iroot != jroot
		root = root_union!(ds, iroot, jroot)
		comp_sizes[root] = comp_sizes[iroot] + comp_sizes[jroot]
		return root
	else
		return iroot
	end
end

function recover_add_nodes(layers, attack_nodes)
	n = nv(layers[1])
	max_comp_sizes = Int[]
	ds = IntDisjointSets(n)
	comp_sizes = ones(Int, n)
	max_comp_size = 1
	presents = [fill(true, n) for _ in 1:length(layers)]
	for attack_node in attack_nodes
		presents[attack_node.layer][attack_node.id] = false
	end

	# add edges between present nodes and record component sizes
	for (layer, present) in zip(layers, presents)
		for e in edges(layer)
			i, j = src(e), dst(e)
			if present[i] && present[j]
				newroot = merge_nodes!(ds, comp_sizes, i, j)
				if comp_sizes[newroot] > max_comp_size
					max_comp_size = comp_sizes[newroot]
				end
			end
		end
	end

	for attack_node in reverse(attack_nodes)
		for j in neighbors(layers[attack_node.layer], attack_node.id)
			if presents[attack_node.layer][j]
				newroot = merge_nodes!(ds, comp_sizes, attack_node.id, j)
				if comp_sizes[newroot] > max_comp_size
					max_comp_size = comp_sizes[newroot]
				end
			end
		end
		presents[attack_node.layer][attack_node.id] = true
		push!(max_comp_sizes, max_comp_size)
	end
	reverse(max_comp_sizes)
end

function is_cyclic_util(g, v, visited, parent)
	visited[v] = true
	for i in neighbors(g,v)
		if !visited[i]
			is_cyclic_util(g,i,visited,v) && return true
		# If an adjacent vertex is visited and not parent of current vertex,
		# then there is a cycle
		elseif parent != i
			return true
		end
	end
	return false
end

function iscyclic(g)
	visited = fill(false,nv(g))
	for i in vertices(g)
		if !visited[i]
			is_cyclic_util(g,i,visited,-1) && return true
		end
	end
	return false
end

function relabel(g, r)
	index = shuffle(vertices(g))
	num = round(Int, r*nv(g))
	if num < 1
		node_map = Dict(zip(index,index))
	elseif num >= nv(g)
		index1 = shuffle(index)
		node_map = Dict(zip(index,index1))
	else
		subindex = shuffle(index[1:num])
		index1 = vcat(subindex, index[num+1:end])
		node_map = Dict(zip(index,index1))
	end
	h = SimpleGraph(nv(g))
	for edge in edges(g)
		add_edge!(h, node_map[src(edge)], node_map[dst(edge)])
	end
	return [g;h]
end

function correlated_multigraph(g1,g2;correlation="UC")
	deg1 = degree(g1)
	deg2 = degree(g2)
	n = nv(g2)
	if correlation == "MP"
		node_map = Dict(zip(sortperm(deg2),sortperm(deg1)))
	elseif correlation == "MN"
		node_map = Dict(zip(sortperm(deg2,rev=true),sortperm(deg1)))
	else
		node_map = Dict(zip(shuffle(1:n),shuffle(1:n)))
	end
	g3 = SimpleGraph(nv(g2))
	for edge in edges(g2)
		add_edge!(g3, node_map[src(edge)], node_map[dst(edge)])
	end
	return [g1;g3]
end

function attack2present(layers, attack_nodes)
	presents = [fill(true, nv(layers[1])) for _ in eachindex(layers)]
	for i in attack_nodes
		presents[i.layer][i.id] = false
	end
	return presents
end

function mymean(x,y)
    lx = length(x)
    ly = length(y)
    z = zeros(max(lx,ly))
    if lx < ly
        for i in eachindex(z)
            if i <= lx
                z[i] = (x[i]+y[i])/2
            else
                z[i] = y[i]
            end
        end
        return z
    else
        for i in eachindex(z)
            if i <= ly
                z[i] = (x[i]+y[i])/2
            else
                z[i] = x[i]
            end
        end
        return z
    end
end