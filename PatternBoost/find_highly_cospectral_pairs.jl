#Skips optimization in the case that the score is already positive to avoid the machine "over-optimizing".

#potential optimisations:
    #input form
    #input with delimiters
    #change various randomness/deterministic features
    #change "cospectralise" routine


include("constants.jl")
include("reconstruction_tests.jl")
using JSON
using Polynomials
using DataStructures
using Random
using DeferredAcceptance
using NautyGraphs, Graphs#, LightGraphs
const N = 10 #number of vertices
const p = 0.4 #initialising Erdos-Renyi graph parameter
const k = 1 #level of cospectrality we seek to create

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
            adj[i,j] = line[place]
            adj[j,i] = line[place]
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


function loop_diff(G, H)
    total = 0
    G_tup = spec_tuple(G)
    H_tup = spec_tuple(H)
    for i in eachindex(G_tup)
        total += (G_tup[i]-H_tup[i])^2
    end
    return total
end

function r_set_matching_spec(G, H, r)
    tuples = []
    G_stups = []
    H_stups = []
    for tup in combinations(1:size(G,1), r)
        M_G = G[1:end .∉ [tup], 1:end .∉ [tup]]
        M_H = H[1:end .∉ [tup], 1:end .∉ [tup]]
        push!(G_stups, Any[tup, spec_tuple(M_G)])
        push!(H_stups, Any[tup, spec_tuple(M_H)])
        push!(tuples, tup)
    end

    for gstup in G_stups
        ranking = []
        for hstup in H_stups
            total = 0
            for i in eachindex(gstup[2])
                total += (hstup[2][i]-gstup[2][i])^2
            end
            push!(ranking, [total, hstup[1]])
        end
        sort!(ranking)
        push!(gstup, ranking)
    end

    for gstup in H_stups
        ranking = []
        for hstup in G_stups
            total = 0
            for i in eachindex(gstup[2])
                total += (hstup[2][i]-gstup[2][i])^2
            end
            push!(ranking, [total, hstup[1]])
        end
        sort!(ranking)
        push!(gstup, ranking)
    end

    len = Int(binomial(size(G,1),r))

    students = zeros(Int64, len, len)
    for j in 1:len
        value_seq = G_stups[j][3]
        order = []
        for t in 1:len
            push!(order, value_seq[t][2])
        end
        for i in 1:len
            students[i,j] = findfirst(isequal(tuples[i]), order)
        end
    end

    schools = zeros(Int64, len, len)
    for j in 1:len
        value_seq = H_stups[j][3]
        order = []
        for t in 1:len
            push!(order, value_seq[t][2])
        end
        for i in 1:len
            schools[i,j] = findfirst(isequal(tuples[i]), order)
        end
    end
    capacities = [1 for i in 1:len]
    # schools_tiebroken = STB(schools)
    assignment = DA(students, schools, capacities)[1]
    return assignment
end

function cospectralise(G, H, vertex_selection = "random", edge_addition = "random", edge_deletion = "random")
    #match vertex sets "optimally"
    matching = r_set_matching_spec(G,H,1)

    #pick G or H w.p. 1/2
    setting = rand(1:2)
    if setting == 2
        tmp = copy(G)
        G = copy(H)
        H = tmp
        matching = inverse_permutation(matching)
    end

    len = size(G,1)

    if vertex_selection == "random"
        place = rand(1:len)
    elseif vertex_selection == "optimal"
        #pick the vertex, v, with the maximal loop difference in matching
        max_dif, place = 0, 1
        for i in 1:len
            dif = loop_diff(G[1:end .∉ [i], 1:end .∉ [i]],
                H[1:end .∉ [matching[i]], 1:end .∉ matching[[i]]])
            if dif > max_dif
                max_dif, place = dif, i
            end
        end
    else
        error("edge deletion must be 'random' or 'optimal'")
    end

    #pick and edge adj to v delete
    if edge_deletion == "random"
        deg_v = sum(G[place,:])
        if deg_v == 0
            tmp = nothing
        else
            choice = rand(1:deg_v)
            indices = findall(isequal(1), G[place,:])
            G[indices[choice], place] = 0
            G[place, indices[choice]] = 0
        end
    elseif edge_deletion == "optimal"
        println("to do")
    else
        error("edge deletion must be 'random' or 'optimal'")
    end


    #add edge back randomly/optimally
    if edge_addition == "random"
        anti_deg_v = len - sum(G[place,:])
        if anti_deg_v == 0
            choice = 1
        else
            choice = rand(1:(anti_deg_v-1))
        end
        indices = findall(isequal(0), G[place,:])
        filter!(e -> e != place, indices)
        G[indices[choice], place] = 1
        G[place, indices[choice]] = 1
    elseif edge_addition == "optimal"
        println("to do")
    else
        error("edge addition must be 'random' or 'optimal'")
    end

    #if G and H wer swapped, swap back now
    if setting == 2
        tmp = copy(G)
        G = copy(H)
        H = tmp
    end

    return (G,H)
end

# G = [0 0 1 1 0 1; 0 0 1 0 0 0; 1 1 0 1 0 0; 1 0 1 0 1 0; 0 0 0 1 0 0; 1 0 0 0 0 0]
# H = [0 1 0 1 0 1; 1 0 1 0 0 1; 0 1 0 1 0 0; 1 0 1 0 1 0; 0 0 0 1 0 0; 1 1 0 0 0 0]

# G=adjacency_matrix(erdos_renyi(8,0.6))
# H=adjacency_matrix(erdos_renyi(8,0.6))

# loop_diff(G,H)

# for i in 1:1000
#     G_i,H_i = cospectralise(G,H)
#     if loop_diff(G_i,H_i) < loop_diff(G,H)
#         println(loop_diff(G_i,H_i))
#         # println(loop_diff(G,H))
#         # println(loop_diff(G_i,H_i) <= loop_diff(G,H))
#         G,H=G_i,H_i
#     end
#     if key_cospectrum(G,1)!=key_cospectrum(H,1) && loop_diff(G,H) == 0
#         println(G)
#         println(H)
#         println(loop_diff(G,H))
#         break
#     end
# end

#The local search algorithm for Pattern boost.

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::OBJ_TYPE
    #input a database and a flattened graph pair, return the improved graph similarly

    # Check if current object has positive reward - if so, skip greedy search
    current_reward = reward_calc(obj)
    if current_reward > 0
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
        M_G = G[1:end .∉ [tuples[i]], 1:end .∉ [tuples[i]]]
        M_H = H[1:end .∉ [tuples[matching[i]]], 1:end .∉ [tuples[matching[i]]]]
        score = loop_diff(M_G,M_H)
        push!(ordered_tuples, [score, tuples[i], tuples[matching[i]]])
    end
    sort!(ordered_tuples, rev=true) # potentially better to not reverse

    for i in eachindex(ordered_tuples)
        # M_G = G[1:end .∉ [ordered_tuples[i][2]], 1:end .∉ [ordered_tuples[i][2]]]
        # M_H = H[1:end .∉ [ordered_tuples[i][3]], 1:end .∉ [ordered_tuples[i][3]]]

        rows_MG2 = setdiff(axes(G, 1), ordered_tuples[i][2])
        cols_MG2 = setdiff(axes(G, 2), ordered_tuples[i][2])
        rows_MH2 = setdiff(axes(H, 1), ordered_tuples[i][3])
        cols_MH2 = setdiff(axes(H, 2), ordered_tuples[i][3])

        M_G_2= @view G[rows_MG2, cols_MG2]
        M_H_2= @view G[rows_MH2, cols_MH2]

        M_G_2, M_H_2 = cospectralise(M_G_2, M_H_2)
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
        M_G = G[1:end .∉ [tuples[i]], 1:end .∉ [tuples[i]]]
        M_H = H[1:end .∉ [tuples[matching[i]]], 1:end .∉ [tuples[matching[i]]]]
        score -= loop_diff(M_G,M_H)
    end
    
    return score
end

#Produces a random pair of graphs.
function empty_starting_point()::OBJ_TYPE
    #TODO: recover sample generation   
    G=adjacency_matrix(erdos_renyi(N,p))
    H=adjacency_matrix(erdos_renyi(N,p))
    return flatten_pair(G,H)
end

