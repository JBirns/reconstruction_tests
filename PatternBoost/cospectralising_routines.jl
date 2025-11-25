using JSON
using Polynomials
using DataStructures
using Random
using DeferredAcceptance
using NautyGraphs, Graphs#, LightGraphs

include("graph_line_formatting.jl")
include("reconstruction_tests.jl")

function r_set_matching_spec(G, H, r)
    tuples = []
    G_stups = []
    H_stups = []
    for tup in combinations(1:size(G,1), r)
        rows_MG = collect(setdiff(axes(G, 1), tup))
        cols_MG = collect(setdiff(axes(G, 2), tup))
        rows_MH = collect(setdiff(axes(H, 1), tup))
        cols_MH = collect(setdiff(axes(H, 2), tup))
        M_G = G[rows_MG, cols_MG]
        M_H = H[rows_MH, cols_MH]
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


# function cospectralise(G, H, vertex_selection = "random", edge_addition = "random", edge_deletion = "random")
#     #match vertex sets "optimally"
#     matching = r_set_matching_spec(G,H,1)

#     #pick G or H w.p. 1/2
#     setting = rand(1:2)
#     if setting == 2
#         tmp = copy(G)
#         G = copy(H)
#         H = tmp
#         matching = inverse_permutation(matching)
#     end

#     len = size(G,1)

#     if vertex_selection == "random"
#         place = rand(1:len)
#     elseif vertex_selection == "optimal"
#         #pick the vertex, v, with the maximal loop difference in matching
#         max_dif, place = 0, 1
#         for i in 1:len
#             rows_G = collect(setdiff(axes(G, 1), [i]))
#             cols_G = collect(setdiff(axes(G, 2), [i]))
#             rows_H = collect(setdiff(axes(H, 1), [matching[i]]))
#             cols_H = collect(setdiff(axes(H, 2), [matching[i]]))
#             dif = loop_diff(G[rows_G, cols_G], H[rows_H, cols_H])
#             if dif > max_dif
#                 max_dif, place = dif, i
#             end
#         end
#     else
#         error("edge deletion must be 'random' or 'optimal'")
#     end

#     #pick and edge adj to v delete
#     if edge_deletion == "random"
#         deg_v = sum(G[place,:])
#         if deg_v == 0
#             tmp = nothing
#         else
#             choice = rand(1:deg_v)
#             indices = findall(isequal(1), G[place,:])
#             G[indices[choice], place] = 0
#             G[place, indices[choice]] = 0
#         end
#     elseif edge_deletion == "optimal"
#         println("to do")
#     else
#         error("edge deletion must be 'random' or 'optimal'")
#     end


#     #add edge back randomly/optimally
#     if edge_addition == "random"
#         anti_deg_v = len - sum(G[place,:])
#         indices = findall(isequal(0), G[place,:])
#         filter!(e -> e != place, indices)
#         if isempty(indices)
#             # No edges to add, skip
#         elseif length(indices) == 1
#             choice = 1
#             G[indices[choice], place] = 1
#             G[place, indices[choice]] = 1
#         else
#             choice = rand(1:length(indices))
#             G[indices[choice], place] = 1
#             G[place, indices[choice]] = 1
#         end
#     elseif edge_addition == "optimal"
#         println("to do")
#     else
#         error("edge addition must be 'random' or 'optimal'")
#     end

#     #if G and H were swapped, swap back now
#     if setting == 2
#         tmp = copy(G)
#         G = copy(H)
#         H = tmp
#     end

#     return (G,H)
# end


function cospectralise(G, H, num_changes = 10)
    #pick G or H w.p. 1/2
    setting = rand(1:2)
    if setting == 2
        tmp = copy(G)
        G = copy(H)
        H = tmp
    end

    n = size(G,1)
    G_line = flatten_adj(G)
    m = sum(G_line)

    G_line[randperm(length(G_line))[1:num_changes]] .⊻= 1






    #if G and H were swapped, swap back now
    if setting == 2
        tmp = copy(G)
        G = copy(H)
        H = tmp
    end

    return (G,H)
end


G = [0 0 1 1 0 1; 0 0 1 0 0 0; 1 1 0 1 0 0; 1 0 1 0 1 0; 0 0 0 1 0 0; 1 0 0 0 0 0]
H = [0 1 0 1 0 1; 1 0 1 0 0 1; 0 1 0 1 0 0; 1 0 1 0 1 0; 0 0 0 1 0 0; 1 1 0 0 0 0]

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