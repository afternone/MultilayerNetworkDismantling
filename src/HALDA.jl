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

function HALDA(layers, T=sum(nv.(layers)))
    g = foldl(union, layers)
    n = nv(g)
    mask = fill(false, n)
    neis = zeros(Int, n)
	deg = degree(g)
    attacked_nodes = @NamedTuple{layer::Int,id::Int}[]
	presents = [fill(true, n) for _ in eachindex(layers)]
	degs = degree.(layers)
	H = SortedDict{Tuple{Int,Int},Set{eltype(attacked_nodes)}}()

	for l in eachindex(layers)
		for u in vertices(layers[l])
			degu = (deg[u],degs[l][u])
			if haskey(H, degu)
				push!(H[degu], (layer=l,id=u))
			else
				push!(H, degu=>Set([(layer=l,id=u)]))
			end
		end
	end

	degmax = last(H)[1]
	done = false
	for t in 1:T
		i = rand(H[degmax])
		delete!(H[degmax],i)
		presents[i.layer][i.id] = false
		push!(attacked_nodes, i)

		newdegi = union_degree(layers, presents, i.id, mask, neis)
		if newdegi < deg[i.id]
			for l in eachindex(layers)
				if presents[l][i.id]
					degi = (deg[i.id],degs[l][i.id])
					delete!(H[degi],(layer=l,id=i.id))
					degi = (newdegi,degs[l][i.id])
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
            if presents[i.layer][j]
                newdegj = union_degree(layers, presents, j, mask, neis)
				degj = (deg[j],degs[i.layer][j])
				delete!(H[degj],(layer=i.layer,id=j))
				degs[i.layer][j] -= 1
				degj = (newdegj,degs[i.layer][j])
				if haskey(H, degj)
					push!(H[degj], (layer=i.layer,id=j))
				else
					push!(H, degj=>Set([(layer=i.layer,id=j)]))
				end

				for l in eachindex(layers)
					if l != i.layer && presents[l][j] && newdegj < deg[j]
						degj = (deg[j],degs[l][j])
						delete!(H[degj],(layer=l,id=j))
						degj = (newdegj,degs[l][j])
						if haskey(H, degj)
							push!(H[degj], (layer=l,id=j))
						else
							push!(H, degj=>Set([(layer=l,id=j)]))
						end
					end
                end
                deg[j] = newdegj
			end
		end

		while isempty(H[degmax])
			delete!(H,degmax)
			degmax = last(H)[1]
			if degmax < (1,0)
				done = true
				break
			end
		end
		done && break
	end
	attacked_nodes
end