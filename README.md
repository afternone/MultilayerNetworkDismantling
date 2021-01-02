# MultilayerNetworkDismantling.jl
Find the minimum number of nodes such that after their removal the union graph of multilayer network is broken into connected components of sub-extensive size.

1. You should install Julia first, and then install necessary packages from the Julia REPL as follows.
```
] add LightGraphs DelimitedFiles Random
```
2. Clone this repository.
```
git clone https://github.com/afternone/MultilayerNetworkDismantling.git
```
3. Prepare a multiplex network file (e.g. `toy_multi_net.txt`) with each line represents a layer edge as follows:
```
layer, source, target
1, 1, 2 # an edge between nodes 1 and 2 in layer 1
2, 1, 3 # an edge between nodes 1 and 3 in layer 2
1, 2, 4
1, 1, 5
2, 3, 4
3, 2, 5
```
4. Change to the `src` directory and open terminal from there (don't forget to add Julia to the environment variables `PATH`), and then using the following command to run CoreHLDA on the multiplex file.
```
julia main.jl toy_multi_net.txt
```
4. In the same directory of `toy_multi_net.txt`, you will get the dismantling set file `toy_multi_net.dismantling_set` with each line represents an attacked layer node in format `layer, id`.
