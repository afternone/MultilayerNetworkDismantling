"""
Optimally break tree until the giant component size is less than threshold.
"""
function tree_break(layers, threshold, attacked_nodes=@NamedTuple{layer::Int,id::Int}[])
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
    return attack_nodes
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
