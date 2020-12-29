function CoreHLDA(layers, threshold=round(Int,sqrt(nv(layers[1]))))
	decore_nodes = decore(layers)
	treebreak_nodes = tree_break(layers, threshold, decore_nodes)
	attack_nodes = vcat(decore_nodes, treebreak_nodes)
	reinsert_nodes, presents = reverse_greedy!(layers, threshold, attack_nodes)
	#all_attack_nodes = [i for i in [decore_nodes;treebreak_nodes] if !presents[i.layer][i.id]]
	return setdiff(attack_nodes, reinsert_nodes)
end