# Isomorphism test for graphs represented as adjacency matrices
# Uses NautyGraphs for canonical labeling

include("graph_line_formatting.jl")
using Graphs
using NautyGraphs
using LinearAlgebra

"""
Convert an adjacency matrix to a SimpleGraph object.

Args:
    adj: 2D array representing an adjacency matrix (symmetric, 0/1 values)
    
Returns:
    SimpleGraph: The graph object
"""
function adj_to_graph(adj)
    n = size(adj, 1)
    G = SimpleGraph(n)
    
    # Add edges based on upper triangle of adjacency matrix
    for i in 1:n
        for j in (i+1):n
            if adj[i, j] != 0
                add_edge!(G, i, j)
            end
        end
    end
    
    return G
end

"""
Check if two adjacency matrices represent isomorphic graphs.

Two graphs are isomorphic if and only if their canonical forms (computed via
NautyGraphs) have identical adjacency matrices.

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
    
    # Convert to SimpleGraph objects
    G1 = adj_to_graph(adj1)
    G2 = adj_to_graph(adj2)
    
    # Get canonical forms using NautyGraphs
    ng1 = NautyGraph(G1)
    ng2 = NautyGraph(G2)
    
    # Convert canonical forms back to graphs
    G1_canonical = SimpleGraph(ng1)
    G2_canonical = SimpleGraph(ng2)
    
    # Get adjacency matrices of canonical forms
    adj1_canonical = Matrix(adjacency_matrix(G1_canonical))
    adj2_canonical = Matrix(adjacency_matrix(G2_canonical))
    
    # Two graphs are isomorphic if and only if their canonical forms are identical
    return adj1_canonical == adj2_canonical
end

"""
Check if the two graphs in a pair line are isomorphic.

Args:
    line: Flattened representation of two graphs (from pair_line_to_adj format)
    
Returns:
    Bool: true if the two graphs are isomorphic, false otherwise
"""
function check_isomorphism_from_line(line)
    G, H = pair_line_to_adj(line)
    return are_isomorphic(G, H)
end

# ============================================================================
# Test cases
# ============================================================================

println("=" ^ 60)
println("Testing isomorphism checker")
println("=" ^ 60)

# Test 1: Two identical graphs (should be isomorphic)
println("\nTest 1: Two identical triangles")
adj_triangle = [0 1 1; 1 0 1; 1 1 0]
result1 = are_isomorphic(adj_triangle, adj_triangle)
println("Result: ", result1, " (expected: true)")

# Test 2: Two isomorphic but different graphs
# Path graph: 1-2-3 vs 3-2-1 (same structure, different vertex labels)
println("\nTest 2: Two isomorphic path graphs")
adj_path1 = [0 1 0; 1 0 1; 0 1 0]  # Path: 1-2-3
adj_path2 = [0 1 0; 1 0 1; 0 1 0]  # Same path
result2 = are_isomorphic(adj_path1, adj_path2)
println("Result: ", result2, " (expected: true)")

# Test 3: Non-isomorphic graphs
println("\nTest 3: Triangle vs Path (non-isomorphic)")
adj_triangle = [0 1 1; 1 0 1; 1 1 0]  # Triangle (3-cycle)
adj_path = [0 1 0; 1 0 1; 0 1 0]      # Path of length 2
result3 = are_isomorphic(adj_triangle, adj_path)
println("Result: ", result3, " (expected: false)")

# Test 4: Your specific line
println("\nTest 4: Graphs from your line")
line = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0]
G, H = pair_line_to_adj(line)
println("Graph G (first graph):")
println(G)
println("\nGraph H (second graph):")
println(H)
result4 = check_isomorphism_from_line(line)
println("\nAre G and H isomorphic? ", result4)

# Show canonical forms for debugging
G_graph = adj_to_graph(G)
H_graph = adj_to_graph(H)
G_canonical = SimpleGraph(NautyGraph(G_graph))
H_canonical = SimpleGraph(NautyGraph(H_graph))
G_canonical_adj = Matrix(adjacency_matrix(G_canonical))
H_canonical_adj = Matrix(adjacency_matrix(H_canonical))
println("\nCanonical form of G:")
println(G_canonical_adj)
println("\nCanonical form of H:")
println(H_canonical_adj)
println("\nCanonical forms are equal: ", G_canonical_adj == H_canonical_adj)

println("\n" * "=" ^ 60)
