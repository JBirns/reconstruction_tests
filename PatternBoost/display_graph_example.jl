# Example: How to display a graph in Julia from an adjacency matrix

# Method 1: Using Graphs.jl and GraphPlot.jl (most common)
include("graph_line_formatting.jl")
using Graphs
using GraphPlot

# Create a graph from an adjacency matrix
function display_graph_from_adj(adj_matrix)
    # Convert adjacency matrix to a Graphs.jl graph
    n = size(adj_matrix, 1)
    G = SimpleGraph(n)
    
    # Add edges based on adjacency matrix
    for i in 1:n
        for j in (i+1):n
            if adj_matrix[i, j] != 0
                add_edge!(G, i, j)
            end
        end
    end
    
    # Plot the graph
    # Layout options: spring_layout, circular_layout, tree_layout, shell_layout, etc.
    gplot(G, 
          layout=spring_layout,  # or circular_layout, tree_layout, shell_layout, etc.
          nodefillc="lightblue",
          nodesize=0.1,
          edgestrokec="red",
          nodelabel=1:nv(G))
end

line = [1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1]   
display_graph_from_adj(pair_line_to_adj(line)[1])
display_graph_from_adj(pair_line_to_adj(line)[2])

