function HMD(layers)
    attack_nodes = @NamedTuple{layer::Int,id::Int}[]
    g = foldl(union, layers)
	deg = degree(g)
	degmax = maximum(deg)
	H = [Set{Int}() for i=1:degmax]

	for i in vertices(g)
		deg[i] > 0 && push!(H[deg[i]], i)
	end

	while degmax > 0
		if !isempty(H[degmax])
			i = rand(H[degmax])
			delete!(H[degmax], i)
            deg[i] = 0
            for l in eachindex(layers)
                push!(attack_nodes, (layer=l,id=i))
            end
		else
			degmax -= 1
		end
	end
	return attack_nodes
end
