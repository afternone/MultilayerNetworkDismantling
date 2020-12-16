using LightGraphs
using DataStructures
using Base.Iterators

function union_degree(layers, presents, i)
	nei = Set{Int}()
	for l in eachindex(layers)
		if presents[l][i]
			for j in neighbors(layers[l],i)
				if presents[l][j]
					push!(nei,j)
				end
			end
		end
	end
	return length(nei)
end

function merge_layer(layers, presents)
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

function decored_graph(layers, attack_nodes)
	presents = [fill(true, nv(layers[1])) for _ in eachindex(layers)]
	for i in attack_nodes
		presents[i.layer][i.id] = false
	end
	merge_layer(layers, presents)
end

function MHSDA(layers, attacked_nodes=@NamedTuple{layer::Int,id::Int}[])
	attack_nodes = @NamedTuple{layer::Int,id::Int}[]
	presents = [fill(true, nv(layers[1])) for _ in 1:length(layers)]
	for attack_node in attacked_nodes
		presents[attack_node.layer][attack_node.id] = false
	end
	g = merge_layer(layers, presents)
	merge_degs = degree(g)

	# initialize degrees
	layer_degs = degree.(layers)
	for i in eachindex(layers)
		g, present = layers[i], presents[i]
		for u in vertices(g)
			if present[u]
				for v in neighbors(g, u)
					if !present[v]
						layer_degs[i][u] -= 1
					end
				end
			else
				layer_degs[i][u] = -1
			end
		end
	end

	degmax = maximum(flatten(map(x->x.+merge_degs,layer_degs)))
	H = [Set{eltype(attack_nodes)}() for i=1:degmax]
	for i in eachindex(layers)
		for u in vertices(layers[i])
			if presents[i][u] && layer_degs[i][u] > 0
				push!(H[merge_degs[u]+layer_degs[i][u]], (layer=i,id=u))
			end
		end
	end

	while degmax > 0
		if !isempty(H[degmax])
			i = rand(H[degmax])
			delete!(H[degmax], i)
			presents[i.layer][i.id] = false
			layer_degs[i.layer][i.id] = -1
			push!(attack_nodes, i)
		else
			degmax -= 1
			continue
		end
		# update neighbors
		for v in neighbors(layers[i.layer], i.id)
			if presents[i.layer][v] && layer_degs[i.layer][v] > 0
				delete!(H[merge_degs[v]+layer_degs[i.layer][v]], (layer=i.layer,id=v))
				layer_degs[i.layer][v] -= 1
				nei = Set{Int}()
				for l in eachindex(layers)
					for w in neighbors(layers[l], v)
						if presents[l][w] && layer_degs[l][w] > 0
							push!(nei, w)
						end
					end
				end
				merge_degs[v] = length(nei)
				layer_degs[i.layer][v] > 0 && push!(H[merge_degs[v]+layer_degs[i.layer][v]], (layer=i.layer,id=v))
			end
		end
	end
	return attack_nodes
end

function MHDA(layers, attacked_nodes=@NamedTuple{layer::Int,id::Int}[])
	attack_nodes = @NamedTuple{layer::Int,id::Int}[]
	presents = [fill(true, nv(layers[1])) for _ in 1:length(layers)]
	for attack_node in attacked_nodes
		presents[attack_node.layer][attack_node.id] = false
	end

	# initialize degrees
	degs = degree.(layers)
	for i in eachindex(layers)
		g, present = layers[i], presents[i]
		for u in vertices(g)
			if present[u]
				for v in neighbors(g, u)
					if !present[v]
						degs[i][u] -= 1
					end
				end
			else
				degs[i][u] = -1
			end
		end
	end

	degmax = maximum(flatten(degs))
	H = [Set{eltype(attack_nodes)}() for i=1:degmax]
	for i in eachindex(layers)
		for u in vertices(layers[i])
			if presents[i][u] && degs[i][u] > 0
				push!(H[degs[i][u]], (layer=i,id=u))
			end
		end
	end

	while degmax > 0
		if !isempty(H[degmax])
			i = rand(H[degmax])
			delete!(H[degmax], i)
			degs[i.layer][i.id] = -1
			push!(attack_nodes, i)
		else
			degmax -= 1
			continue
		end
		# update neighbors
		for v in neighbors(layers[i.layer], i.id)
			if degs[i.layer][v] > 0
				delete!(H[degs[i.layer][v]], (layer=i.layer,id=v))
				degs[i.layer][v] -= 1
				degs[i.layer][v] > 0 && push!(H[degs[i.layer][v]], (layer=i.layer,id=v))
			end
		end
	end
	return attack_nodes
end

function multilayer_decore(layers, k=2)
    g = foldl(union, layers)
    n = nv(g)
	deg = degree(g)
    attacked_nodes = @NamedTuple{layer::Int,id::Int}[]
	presents = [fill(true, n) for _ in eachindex(layers)]
	degs = degree.(layers)
	H = SortedDict{Tuple{Int,Int},Set{eltype(attacked_nodes)}}()
	push!(H, (k-1,k-1)=>Set{eltype(attacked_nodes)}())
	for l in eachindex(layers)
		for u in vertices(layers[l])
			degu = max((deg[u],degs[l][u]),(k-1,k-1))
			if haskey(H, degu)
				push!(H[degu], (layer=l,id=u))
			else
				push!(H, degu=>Set([(layer=l,id=u)]))
			end
		end
	end
	degmax = last(H)[1]
	d = isempty(H[(k-1,k-1)]) ? degmax : (k-1,k-1)
	for t in 1:length(layers)*n
		i = rand(H[d])
		delete!(H[d],i)
		presents[i.layer][i.id] = false
		d > (k-1,k-1) && push!(attacked_nodes, i)

		newdegi = union_degree(layers, presents, i.id)
		if d > (k-1,k-1) && newdegi < deg[i.id]
			for l in eachindex(layers)
				if presents[l][i.id]
					degi = max((deg[i.id],degs[l][i.id]), (k-1,k-1))
					delete!(H[degi],(layer=l,id=i.id))
					deg[i.id] = newdegi
					degi = max((deg[i.id],degs[l][i.id]), (k-1,k-1))
					if haskey(H, degi)
						push!(H[degi], (layer=l,id=i.id))
					else
						push!(H, degi=>Set([(layer=l,id=i.id)]))
					end
				end
			end
		end

		for j in neighbors(layers[i.layer], i.id)
			newdegj = union_degree(layers, presents, j)

			if presents[i.layer][j] && deg[j] ≥ k
				degj = max((deg[j],degs[i.layer][j]),(k-1,k-1))
				delete!(H[degj],(layer=i.layer,id=j))
				degs[i.layer][j] -= 1
				degj = max((newdegj,degs[i.layer][j]),(k-1,k-1))
				if haskey(H, degj)
					push!(H[degj], (layer=i.layer,id=j))
				else
					push!(H, degj=>Set([(layer=i.layer,id=j)]))
				end

				for l in eachindex(layers)
					if l != i.layer && presents[l][j] && newdegj < deg[j]
						degj = max((deg[j],degs[l][j]),(k-1,k-1))
						delete!(H[degj],(layer=l,id=j))
						degj = max((newdegj,degs[l][j]),(k-1,k-1))
						if haskey(H, degj)
							push!(H[degj], (layer=l,id=j))
						else
							push!(H, degj=>Set([(layer=l,id=j)]))
						end
					end
				end
			end
			deg[j] = newdegj
		end

		while isempty(H[degmax])
			delete!(H,degmax)
			degmax = last(H)[1]
		end
		d = isempty(H[(k-1,k-1)]) ? degmax : (k-1,k-1)
		degmax < (k,0) && break
	end
	attacked_nodes
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

function ncsize!(S, present, g, i, j)
	!present[i] && return 0
	S[i] != 0 && error("the graph is not acyclic")
	S[i] = 1
	for k in neighbors(g,i)
		if k != j && present[k]
			S[i] += ncsize!(S, present, g, k, i)
		end
	end
	return S[i]
end

"""
Optimally break tree until the giant component size is less than threshold.
"""
function treebreak(layers, threshold, attacked_nodes=@NamedTuple{layer::Int,id::Int}[])
	presents = [fill(true, nv(layers[1])) for _ in 1:length(layers)]
	for attack_node in attacked_nodes
		presents[attack_node.layer][attack_node.id] = false
	end

	# initialize degrees
	degs = degree.(layers)
	for l in eachindex(layers)
		for u in vertices(layers[l])
			if presents[l][u]
				for v in neighbors(layers[l], u)
					if !presents[l][v]
						degs[l][u] -= 1
					end
				end
			else
				degs[l][u] = -1
			end
		end
	end

	g = merge_layer(layers, presents)
	present = foldl((x1,x2)->x1 .| x2, presents)

    S = zeros(Int, nv(g))
    H = MutableBinaryMaxHeap{Tuple{Int64,Int64}}()
	for u in vertices(g)
		if present[u] && S[u] == 0
			push!(H, (ncsize!(S, present, g, u, nothing),u))
		end
	end

	attack_nodes = @NamedTuple{layer::Int,id::Int}[]
    scomps = Int[]
    while !isempty(H)
        scomp, i = pop!(H)
        sender = nothing
        while true
            sizes = [(S[k],k) for k in neighbors(g,i) if k != sender && present[k]]
            isempty(sizes) && break
            M, largest = maximum(sizes)
            if M <= scomp/2
				num, maxl, maxdeg = 0, 0, 0
				for l in eachindex(layers)
					if presents[l][i]
						num += 1
						if degs[l][i] ≥ maxdeg
							maxdeg = degs[l][i]
							maxl = l
						end
					end
				end

				if num > 1
					S[i] = 1
					for k in neighbors(g,i)
						if present[k]
							exist_edge = false
							for l in eachindex(layers)
								if l != maxl && has_edge(layers[l], k, i) && presents[l][k] && presents[l][i]
									exist_edge = true
									break
								end
							end
							if exist_edge
								S[i] += S[k]
							end
						end
					end
					push!(H, (S[i],i))
					for k in neighbors(layers[maxl],i)
						exist_edge = false
						for l in eachindex(layers)
							if l != maxl && has_edge(layers[l], k, i) && presents[l][k] && presents[l][i]
								exist_edge = true
								break
							end
						end
						if presents[maxl][k] && !exist_edge
							if S[k] > 1
								push!(H, (S[k],k))
							end
							rem_edge!(g,i,k)
						end
	                end
				else
					present[i] = false
					for k in neighbors(layers[maxl],i)
	                    if presents[maxl][k]
	                        S[k] > 1 && push!(H, (S[k],k))
							rem_edge!(g,i,k)
	                    end
	                end
				end
                presents[maxl][i] = false
                push!(attack_nodes, (layer=maxl,id=i))
				for k in neighbors(layers[maxl],i)
					if presents[maxl][k]
						degs[maxl][k] -= 1
					end
				end
                push!(scomps, scomp)
                break
            end
            S[i] = 1
            for k in neighbors(g,i)
                if k != largest && present[k]
                    S[i] += S[k]
                end
            end
            sender, i = i, largest
        end
		scomp <= threshold && break
    end
    return attack_nodes, scomps
end

function correlated_multigraph(g1,g2;correlation="MP")
	deg1 = degree(g1)
	deg2 = degree(g2)
	n = nv(g2)
	if correlation == "MP"
		node_map = Dict(zip(sortperm(deg2),sortperm(deg1)))
	elseif correlation == "MN"
		node_map = Dict(zip(sortperm(deg2,rev=true),sortperm(deg1)))
	else
		node_map = Dict(zip(rand(1:n,n),rand(1:n,n)))
	end
	g3 = SimpleGraph(nv(g2))
	for edge in edges(g2)
		add_edge!(g3, node_map[src(edge)], node_map[dst(edge)])
	end
	g3
end


