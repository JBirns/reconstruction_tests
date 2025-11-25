#Skips optimization in the case that the score is already 0 to avoid the machine "over-optimizing".

#potential optimisations:


include("constants.jl")
include("reconstruction_tests.jl")
include("cospectralising_routines.jl")
include("graph_line_formatting.jl")
using JSON
using Polynomials
using DataStructures
using Random
using LinearAlgebra
using DeferredAcceptance
using NautyGraphs, Graphs#, LightGraphs
const N = 10 #number of vertices
const M = 20 #initialising #edges
const k = 1 #level of cospectrality we seek to create
const cospectralising_number = 10 # number of times we iterate cospectralising

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


#The local search algorithm for Pattern boost.

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::OBJ_TYPE
    #input a database and a flattened graph pair, return the improved graph similarly

    # Check if current object has max reward - if so, skip greedy search
    current_reward = reward_calc(obj)
    if current_reward == 0
        return obj  # Return the original object unchanged
    end

    line = obj
    G,H = pair_line_to_adj(line)
    n = size(G,1)
    
    tuples = []
    for tup in combinations(1:n, k)
        push!(tuples, tup)
    end

    matching = r_set_matching_spec(G, H, k)

    ordered_tuples = []
    for i in eachindex(tuples)
        rows_MG = collect(setdiff(axes(G, 1), tuples[i]))
        cols_MG = collect(setdiff(axes(G, 2), tuples[i]))
        rows_MH = collect(setdiff(axes(H, 1), tuples[matching[i]]))
        cols_MH = collect(setdiff(axes(H, 2), tuples[matching[i]]))
        M_G = G[rows_MG, cols_MG]
        M_H = H[rows_MH, cols_MH]
        score = loop_diff(M_G,M_H)
        push!(ordered_tuples, [score, tuples[i], tuples[matching[i]]])
    end
    sort!(ordered_tuples, rev=true) # potentially better to not reverse

    for i in eachindex(ordered_tuples)
        # M_G = G[1:end .∉ [ordered_tuples[i][2]], 1:end .∉ [ordered_tuples[i][2]]]
        # M_H = H[1:end .∉ [ordered_tuples[i][3]], 1:end .∉ [ordered_tuples[i][3]]]

        rows_MG2 = collect(setdiff(axes(G, 1), ordered_tuples[i][2]))
        cols_MG2 = collect(setdiff(axes(G, 2), ordered_tuples[i][2]))
        rows_MH2 = collect(setdiff(axes(H, 1), ordered_tuples[i][3]))
        cols_MH2 = collect(setdiff(axes(H, 2), ordered_tuples[i][3]))

        # Skip if submatrix is empty (all vertices removed)
        if isempty(rows_MG2) || isempty(rows_MH2)
            continue
        end

        # Create copies since cospectralise creates new matrices when swapping
        M_G_2 = copy(G[rows_MG2, cols_MG2])
        M_H_2 = copy(H[rows_MH2, cols_MH2])

        for i in 1:cospectralising_number
            G_i,H_i = cospectralise(M_G_2, M_H_2)
            if loop_diff(G_i,H_i) < loop_diff(M_G_2,M_H_2)
                M_G_2, M_H_2 = G_i,H_i
            end
        end
        
        # Copy results back to original matrices
        for (idx1, r) in enumerate(rows_MG2)
            for (idx2, c) in enumerate(cols_MG2)
                G[r, c] = M_G_2[idx1, idx2]
            end
        end
        for (idx1, r) in enumerate(rows_MH2)
            for (idx2, c) in enumerate(cols_MH2)
                H[r, c] = M_H_2[idx1, idx2]
            end
        end
    end

    return flatten_pair(G,H)
end

#Computes the reward function for the input graph
function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    line = obj
    G, H = pair_line_to_adj(line)
    n = size(G,1)


    if NautyGraph(G) == NautyGraph(H) #####THIS IS NOT A FOOLPROOF ISOMORPHISM TEST
        return -1e9  # punishment for isomorphism
    end

    matching = r_set_matching_spec(G, H, k)
    tuples = []
    score = 0

    for tup in combinations(1:n, k)
        push!(tuples, tup)
    end
    for i in eachindex(tuples)
        rows_MG = collect(setdiff(axes(G, 1), tuples[i]))
        cols_MG = collect(setdiff(axes(G, 2), tuples[i]))
        rows_MH = collect(setdiff(axes(H, 1), tuples[matching[i]]))
        cols_MH = collect(setdiff(axes(H, 2), tuples[matching[i]]))
        M_G = G[rows_MG, cols_MG]
        M_H = H[rows_MH, cols_MH]
        score -= loop_diff(M_G,M_H)
    end
    
    #punish disconnectivity
    if is_connected(adj_to_graph(G))==false || is_connected(adj_to_graph(H))==false
        score -= 1e5    
    end

    return score
end

#Produces a random pair of graphs.
function empty_starting_point()::OBJ_TYPE
    #TODO: recover sample generation   
    G=line_to_adj(shuffle([ones(Int, M); zeros(Int, Int(N*(N-1)/2)-M)]))
    H=line_to_adj(shuffle([ones(Int, M); zeros(Int, Int(N*(N-1)/2)-M)]))
    return flatten_pair(G,H)
end

