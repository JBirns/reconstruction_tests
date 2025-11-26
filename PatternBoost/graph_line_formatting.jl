using LinearAlgebra
include("constants.jl")
using Random

function line_to_adj_no_delimiter(line)
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

function random_pair()
    M = rand(1:Int(N*(N-1)/2))
    G=line_to_adj_no_delimiter(shuffle([ones(Int, M); zeros(Int, Int(N*(N-1)/2)-M)]))
    H=line_to_adj_no_delimiter(shuffle([ones(Int, M); zeros(Int, Int(N*(N-1)/2)-M)]))
    return flatten_pair(G,H)
end

function flatten_adj(G)
    line = []
    n = size(G,1)
    for i in 1:(n-1)
        for j in (i+1):n
            push!(line, G[i,j])
        end
        if i != n-1 #delimiter
            push!(line, 2)
        end
    end
    return line
end

function flatten_pair(G,H)
    line1 = flatten_adj(G)
    push!(line(3))
    line2 = flatten_adj(H)
    return [line1;line2]
end

function line_to_adj(line)
    l = length(line) 
    # n = Int((1+sqrt(1+8*l))/2)  no delimiters
    n = Int((-1+sqrt(17+8*l))/2)

    #check has right delimiters
    delimiter_places = []
    i=0
    place = 0
    while n-i >= 3
        place += n-i
        push!(delimiter_places, place)
        i += 1
    end
    for i in eachindex(line)
        if line[i] != 2 && i in delimiter_places
            return random_pair()
        end
        if line[i] == 2 && !(i in delimiter_places)
            return random_pair
        end
    end

    #remove all delimiters
    filter!(x -> x != 2, line)

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
    for i in [1:Int((l-1)/2);Int((l+3)/2):l]
        if line[i] == 3
            return random_pair()
        end
    end
    if line[Int((l+1)/2)] != 3
        return random_pair()
    end
    try
        half = Int((l-1)/2)
        G_line, H_line = [], []
        for i in 1:half
            push!(G_line, line[i])
            push!(H_line, line[i+half])
        end
        return (line_to_adj(G_line), line_to_adj(H_line))
    catch
        return random_pair()
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