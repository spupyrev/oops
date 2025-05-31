# OOPS
This repository contains the source code of a solver for 1-planarity.

OOPS (Optimized One-Planarity Solver) is a practical heuristic for recognizing 
1-planar graphs and several of their important subclasses. A graph
is [1-planar](https://en.wikipedia.org/wiki/1-planar_graph) if it can be drawn in 
the plane such that each edge is crossed at most once. Although recognizing 1-planarity
is NP-complete, OOPS is capable to compute 1-planar drawings of medium sized graphs
within seconds.

## Requirements
Your machine needs to have a C++ compiler with C++17 support. We use GNU Make to build the executables.

## Building & Running
1. Compile the binary by running:

        make

2. Run the executable:

        ./oops -i=data/test4.dot -o=test4.txt

    The tool accepts graphs in the [DOT](https://en.wikipedia.org/wiki/DOT_(graph_description_language)), [GML](https://en.wikipedia.org/wiki/Graph_Modelling_Language), and [g6/s6](https://users.cecs.anu.edu.au/~bdm/data/formats.html) list formats. The output is either in the 
    textual (a list of planar and crossing edges) or visual form (using GML or SVG formats).

For the list of supported options use:

        ./oops -help
	
3. For large instances, it is possible to use an external (e.g., parallel) SAT solver. Download and setup a SAT solver such as [Lingeling](http://fmv.jku.at/lingeling) or [Painless](https://www.lrde.epita.fr/wiki/Painless). Then create a CNF for 1-planarity: 

        ./oops -i=data/test4.dot -dimacs=test4.dimacs

Use the SAT solver to test embeddability of the graph:          

        painless test4.dimacs > result.dimacs

Finally, print the resulting layout:        

        ./oops -i=graphs/test4.dot -dimacs-result=result.dimacs



License
--------
Code is released under the [MIT License](MIT-LICENSE.txt).
