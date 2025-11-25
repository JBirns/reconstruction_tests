using LinearAlgebra

function flatten_adj(G)
    line = []
    n = size(G,1)
    for i in 1:(n-1)
        for j in (i+1):n
            push!(line, G[i,j])
        end
    end
    return line
end

function flatten_pair(G,H)
    line1 = flatten_adj(G)
    line2 = flatten_adj(H)
    return [line1;line2]
end

function line_to_adj(line)
    l = length(line) 
    n = Int((1+sqrt(1+8*l))/2)  
    adj = zeros(Int, n, n) 
    place = 1
    for i in 1:(n-1)
        for j in (i+1):n
            val = line[place]
            # Ensure adjacency matrix only contains 0 or 1
            adj[i,j] = val != 0 ? 1 : 0
            adj[j,i] = val != 0 ? 1 : 0
            place += 1
        end
    end
    return adj
end

function pair_line_to_adj(line)
    l = length(line)
    try
        half = Int(l/2)
        G_line, H_line = [], []
        for i in 1:half
            push!(G_line, line[i])
            push!(H_line, line[i+half])
        end
        return (line_to_adj(G_line), line_to_adj(H_line))
    catch
        return (adjacency_matrix(erdos_renyi(N,p)),adjacency_matrix(erdos_renyi(N,p)))
    end
end

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

# println(pair_line_to_adj([0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0]))