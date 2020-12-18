"""
High Layer Degree Adaptive
"""
function HLDA(layers, attacked_nodes=@NamedTuple{layer::Int,id::Int}[])
	attack_nodes = @NamedTuple{layer::Int,id::Int}[]
	presents = [fill(true,nv(layers[1])) for _ in eachindex(layers)]
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

	degmax = maximum(Iterators.flatten(degs))
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
			u = rand(H[degmax])
			delete!(H[degmax], u)
			degs[u.layer][u.id] = -1
			push!(attack_nodes, u)
		else
			degmax -= 1
			continue
		end
		# update neighbors
		for v in neighbors(layers[u.layer], u.id)
			if degs[u.layer][v] > 0
				delete!(H[degs[u.layer][v]], (layer=u.layer,id=v))
				degs[u.layer][v] -= 1
				degs[u.layer][v] > 0 && push!(H[degs[u.layer][v]], (layer=u.layer,id=v))
			end
		end
	end
	return attack_nodes
end