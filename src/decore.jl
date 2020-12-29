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

function decore(layers, k=2)
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
	done = false
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
					degi = max((newdegi,degs[l][i.id]), (k-1,k-1))
					if haskey(H, degi)
						push!(H[degi], (layer=l,id=i.id))
					else
						push!(H, degi=>Set([(layer=l,id=i.id)]))
					end
				end
			end
		end
		deg[i.id] = newdegi

		for j in neighbors(layers[i.layer], i.id)
			newdegj = union_degree(layers, presents, j)

			if presents[i.layer][j]
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
			if degmax ≤ (k-1,k-1)
				done = true
				break
			end
		end
		done && break
		d = isempty(H[(k-1,k-1)]) ? degmax : (k-1,k-1)
	end
	attacked_nodes
end