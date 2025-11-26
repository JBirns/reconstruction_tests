#Skips optimization in the case that the score is already 0 to avoid the machine "over-optimizing".

#potential optimisations:


include("constants.jl")
include("reconstruction_tests.jl")
include("cospectralising_routines.jl")
include("graph_line_formatting.jl")
include("check_isomorphism.jl")
using JSON
using Polynomials
using DataStructures
using Random
using LinearAlgebra
using DeferredAcceptance
using Graphs
using SimpleGraphAlgorithms
# const M = 4 #initialising #edges
const k = 0 #level of cospectrality we seek to create
const cospectralising_number = 10 # number of times we iterate cospectralising


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


    # Use SimpleGraphAlgorithms for reliable isomorphism testing
    if are_isomorphic(G, H)
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
    
    # #punish disconnectivity
    # if Graphs.is_connected(adj_to_graph(G))==false || Graphs.is_connected(adj_to_graph(H))==false
    #     score -= 1e5   
    #     score = score * 2 
    # end

    return score
end

#Produces a random pair of graphs.
function empty_starting_point()::OBJ_TYPE
    #TODO: recover sample generation
    G, H = random_pair()   
    return flatten_pair(G, H)
end

