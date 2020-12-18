module MultilayerNetworkDismantling

using LightGraphs
using DataStructures
using Base.Iterators

export
HLD,
HLDA,
HMD,
HMDA,
HLCIA,
HMCIA,
CoreHLDA,
decore,
tree_break,
reverse_greedy!,
recover_add_nodes

include("decore.jl")
include("tree_break.jl")
include("reverse_greedy.jl")
include("CoreHLDA.jl")
include("HLD.jl")
include("HLDA.jl")
include("HMD.jl")
include("HMDA.jl")
include("HLCIA.jl")
include("HMCIA.jl")
include("utils.jl")

end # module