module MultilayerNetworkDismantling

using LightGraphs
using DataStructures
using Base.Iterators
using StatsBase
using Random

export
OAS,
LCC,
EMD,
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
recover_add_nodes,
mymean

include("decore.jl")
include("tree_break.jl")
include("reverse_greedy.jl")
include("CoreHLDA.jl")
include("OAS.jl")
include("EMD.jl")
include("HLD.jl")
include("HLDA.jl")
include("HMD.jl")
include("HMDA.jl")
include("HLCIA.jl")
include("HMCIA.jl")
include("utils.jl")

end # module