function reverse_greedy!(layers, threshold, attack_nodes)
	threshold = max(1, round(Int, threshold))
	n = nv(layers[1])
	ds = IntDisjointSets(n)
	reinsert_nodes = @NamedTuple{layer::Int,id::Int}[]
	comp_sizes = ones(Int, n) # size of the component where the node in
	mask = fill(false, n) # auxiliary vector
	compos = zeros(Int, n) # auxiliary vector
	max_comp_size = 1 # record giant component size
	present = fill(false, n)
	presents = [fill(true, nv(layers[1])) for _ in eachindex(layers)]
	for attack_node in attack_nodes
		presents[attack_node.layer][attack_node.id] = false
	end

    # add edges between present nodes and record component sizes
    for l in eachindex(layers)
        for i in vertices(layers[l])
			if presents[l][i]
				present[i] = true
                for j in neighbors(layers[l],i)
                    if presents[l][j]
                        newroot = merge_nodes!(ds, comp_sizes, i, j)
                        if comp_sizes[newroot] > max_comp_size
                            max_comp_size = comp_sizes[newroot]
                        end
                    end
                end
            end
        end
	end
	
	smin = threshold
	ibest = (layer=0,id=0)
	while smin ≤ threshold
		smin = threshold
		for l in eachindex(layers)
			for i in vertices(layers[l])
				if !presents[l][i]
					s = present[i] ? 0 : 1
					ncomp = 0
					for j in neighbors(layers[l],i)
						if presents[l][j]
							jroot = find_root(ds,j)
							if !mask[jroot]
								mask[jroot] = true
								ncomp += 1
								compos[ncomp] = jroot
								s += comp_sizes[jroot]
							end
						end
					end
					# reset the auxiliary vector mask
					for k in 1:ncomp
						mask[compos[k]] = false
					end
					s > threshold && continue
					if s < smin
						smin = s
						ibest = (layer=l,id=i)
					end
				end
			end
		end
		ibest == (layer=0,id=0) && break

		for j in neighbors(layers[ibest.layer], ibest.id)
			if presents[ibest.layer][j]
				newroot = merge_nodes!(ds, comp_sizes, ibest.id, j)
				if comp_sizes[newroot] > max_comp_size
					max_comp_size = comp_sizes[newroot]
				end
			end
		end
		max_comp_size > threshold && break
		push!(reinsert_nodes, ibest)
		presents[ibest.layer][ibest.id] = true
		present[ibest.id] = true
	end 
	reinsert_nodes
end