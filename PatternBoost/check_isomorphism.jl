# Isomorphism test for graphs represented as adjacency matrices
# Uses SimpleGraphAlgorithms for isomorphism checking

include("graph_line_formatting.jl")
using Graphs
using SimpleGraphAlgorithms
using SimpleGraphs
using LinearAlgebra

"""
Convert Graphs.SimpleGraph to SimpleGraphs.UndirectedGraph for use with SimpleGraphAlgorithms.

Args:
    g: Graphs.SimpleGraph object
    
Returns:
    SimpleGraphs.UndirectedGraph: The graph in SimpleGraphs format
"""
function graphs_to_simplegraphs(g::SimpleGraph)
    # Get adjacency matrix from Graphs.SimpleGraph
    adj = adjacency_matrix(g)
    # Construct UndirectedGraph from adjacency matrix
    sg = UndirectedGraph(adj)
    
    return sg
end

"""
Check if two adjacency matrices represent isomorphic graphs.

Uses SimpleGraphAlgorithms.is_iso for reliable isomorphism testing.

Args:
    adj1: First adjacency matrix (2D array)
    adj2: Second adjacency matrix (2D array)
    
Returns:
    Bool: true if graphs are isomorphic, false otherwise
"""
function are_isomorphic(adj1, adj2)
    # Quick check: must have same number of vertices
    n1 = size(adj1, 1)
    n2 = size(adj2, 1)
    if n1 != n2
        return false
    end
    
    # Convert to Graphs.SimpleGraph objects
    G1 = adj_to_graph(adj1)
    G2 = adj_to_graph(adj2)
    
    # Convert to SimpleGraphs.UndirectedGraph for SimpleGraphAlgorithms
    SG1 = graphs_to_simplegraphs(G1)
    SG2 = graphs_to_simplegraphs(G2)
    
    # Use SimpleGraphAlgorithms to check isomorphism directly
    return is_iso(SG1, SG2)
end

# G,H = pair_line_to_adj([0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0])
# println(are_isomorphic(G,H))