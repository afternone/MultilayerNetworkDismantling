module MultilayerNetworkDismantling

using LightGraphs
using DataStructures
using Base.Iterators

export
HLDA,
CoreHLDA,
decore,
tree_break,
reverse_greedy!,
recover_add_nodes

include("decore.jl")
include("tree_break.jl")
include("reverse_greedy.jl")
include("CoreHLDA.jl")
include("HLDA.jl")
include("utils.jl")

end # module