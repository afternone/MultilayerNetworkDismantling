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