using Combinatorics
using LinearAlgebra

include("graph6_translation.jl")

function inverse_permutation(p::Vector{Int})
    n = length(p)
    p_inv = zeros(Int, n)
    
    for i in 1:n
        p_inv[p[i]] = i
    end
    
    return p_inv
end

function exec(string)
    return eval(Meta.parse(string))
end

"""
    set_in_place!(arr, place::AbstractString, value)

Set arr at the nested index described by place (e.g. "[4][2,1]") to value.
Supports integer tokens and the token "end" (resolved against the current container).
No use of eval; safe for use inside functions.
"""
function set_in_place!(arr, place::AbstractString, value)
    # collect what's inside each [...]
    bracket_contents = [m.captures[1] for m in eachmatch(r"\[(.*?)\]", place)]
    if isempty(bracket_contents)
        throw(ArgumentError("place must contain at least one bracket like \"[3]\""))
    end

    parent = arr
    n = length(bracket_contents)

    # helper: parse one token (integer or "end")
    parse_token(tok::AbstractString, parent_container, dim_index::Int) = begin
        s = strip(tok)
        if occursin(r"^-?\d+$", s)           # integer literal
            return parse(Int, s)
        elseif lowercase(s) == "end"         # 'end' token -> resolve using size
            # If asked for a dimension beyond ndims, size(...,d) returns 1 for arrays;
            # we handle by using size(parent_container, dim_index).
            return size(parent_container, dim_index)
        else
            throw(ArgumentError("unsupported token '$s'. Only integer literals or 'end' are supported."))
        end
    end

    # iterate through bracket groups
    for i in 1:n
        raw = bracket_contents[i]
        parts = split(raw, ",")
        # build tuple of indices for this bracket
        idxs = Tuple(parse_token(parts[j], parent, j) for j in 1:length(parts))
        if i < n
            # step down one level to become parent for next bracket
            parent = getindex(parent, idxs...)   # parent[idxs...]
        else
            # final bracket: perform setindex! on current parent
            setindex!(parent, value, idxs...)
            return arr
        end
    end

    return arr
end

"""
    sort_in_place!(arr, place::AbstractString; kwargs...)

Sort the array-like object located at place (e.g. "[4][2]") in place.
Works inside functions. Supports integer indices and the token "end".
Additional keyword arguments are passed to sort!.
"""
function sort_in_place!(arr, place::AbstractString; kwargs...)
    # Find everything inside [...] groups
    bracket_contents = [m.captures[1] for m in eachmatch(r"\[(.*?)\]", place)]
    if isempty(bracket_contents)
        # throw(ArgumentError("place must contain at least one bracket like \"[3]\""))
        sort!(arr; kwargs...)
        return arr
    end

    parent = arr
    n = length(bracket_contents)

    # helper to parse tokens
    parse_token(tok::AbstractString, parent_container, dim_index::Int) = begin
        s = strip(tok)
        if occursin(r"^-?\d+$", s)
            return parse(Int, s)
        elseif lowercase(s) == "end"
            return size(parent_container, dim_index)
        else
            throw(ArgumentError("unsupported token '$s' — only integers or 'end' are supported."))
        end
    end

    # Walk down to the target sub-array
    for i in 1:n
        raw = bracket_contents[i]
        parts = split(raw, ",")
        idxs = Tuple(parse_token(parts[j], parent, j) for j in 1:length(parts))
        if i < n
            parent = getindex(parent, idxs...)
        else
            # Target element: sort it in place
            sort!(parent[idxs...]; kwargs...)
            return arr
        end
    end

    return arr
end

"""
    get_in_place(arr, place::AbstractString)

Retrieve the value at place (same syntax as set_in_place!).
"""
function get_in_place(arr, place::AbstractString)
    bracket_contents = [m.captures[1] for m in eachmatch(r"\[(.*?)\]", place)]
    if isempty(bracket_contents)
        throw(ArgumentError("place must contain at least one bracket like \"[3]\""))
    end

    parent = arr
    for i in 1:length(bracket_contents)
        raw = bracket_contents[i]
        parts = split(raw, ",")
        idxs = Tuple((strip(p) |> s -> occursin(r"^-?\d+$", s) ? parse(Int, s) :
                      (lowercase(s) == "end" ? size(parent, length(parts)) :
                       throw(ArgumentError("unsupported token '$s'."))))
                     for p in parts)
        parent = getindex(parent, idxs...)
    end
    return parent
end

function spec_tuple(G::AbstractMatrix)
    A = copy(G)
    out = Int[]
    for i in 1:size(G, 1)-1
        A *= G
        push!(out, tr(A))
    end
    return out
end

function empty_tree(n, depth)
    if n < 0 || depth < 0 || !(n == floor(n)) || !(depth == floor(depth))
        throw(ErrorException("Nonnegative integers only"))
    end
    if n == 0 || depth == 0
        return n
    else
        children = []
        for j in 1:n
            push!(children, [empty_tree(n - 1, depth - 1)])
        end
        out = [n, children]
        return out
    end
end

function re_index(lst::Vector{Int})
    # Traverse the list from right to left
    for i in length(lst):-1:1
        # Count how many elements to the left are smaller than lst[i]
        count = 0
        for j in 1:i-1
            if lst[j] < lst[i]
                count += 1
            end
        end
        # Decrease the current term by the count
        lst[i] -= count
    end
    return lst
end

function relabel_tree(tree, vertices::Int, k::Int, bijection)
    for layer in 1:k
        for tup in combinations(1:vertices, layer)  # same as range(vertices) in Python
            # --- Find old colour ---
            altered_tup = re_index(copy(tup))
            tstr = ""
            for p in altered_tup
                tstr = tstr * "[2]" * string([p]) * "[1]"
            end
            if layer < k
                tstr = tstr * "[1]"
            end
            # old_colour = eval(tstr)
            old_colour = get_in_place(tree,tstr)
            new_colour = bijection[old_colour]

            #--- Relabel all permutations of tup ---
            for perm in permutations(tup)
                altered_perm = re_index(copy(perm))
                tstr = ""
                for p in altered_perm
                    tstr = tstr * "[2]" * string([p]) * "[1]"
                end
                if layer < k
                    tstr = tstr * "[1]"
                end
                # tstr = tstr * "= new_colour"
                # exec(tstr)
                set_in_place!(tree, tstr, new_colour)
            end
        end
    end
    return tree
end

function key_cospectrum(G, k::Int, test=spec_tuple)
    if k == 0
        spec = test(G)
        return ([[spec], [1]])
    end

    tree = empty_tree(size(G, 1), k)
    key = [[], []]

    stup = test(G)

    push!(key[1], [stup, 1])
    push!(key[2], stup)
    tree[1] = 1

    for layer in 1:k
        for tup in combinations(1:size(G,1), layer)
            M = G[1:end .∉ [tup], 1:end .∉ [tup]]
            stup = test(M)
            if !(stup in key[2])
                push!(key[1], [stup, length(key[2])+1])
                push!(key[2], stup)
            end

            # Find colour (index)
            colour = findfirst(x -> x == stup, key[2])

            for perm in permutations(tup)
                if isa(perm, Int) #convert tuples of length 1 to lists
                    perm = [perm]
                end
                tstr = ""
                altered_perm = re_index(copy(perm))
                for p in altered_perm
                    tstr = tstr * "[2]" * string([p]) * "[1]"
                end
                if layer < k #index???
                    tstr = tstr * "[1]"
                end
                set_in_place!(tree, tstr, colour) #used to be exec(tstr * "= colour"), and tstr began with "tree"
            end
        end
    end

    # --- Canonical relabeling ---
    tmp = deepcopy(key[2])
    sort!(tmp; by = x -> (-length(x),x))  #sort lex then by reverse length
    key[2] = tmp

    bijection = Int[]
    for i in 1:length(key[2]) #2:...???? instead
        pos = findfirst(x -> x[1] == key[2][i], key[1])#???
        push!(bijection, key[1][pos][2])
    end

    bijection_perm = inverse_permutation(bijection)
    tree = relabel_tree(tree, size(G,1), k, bijection_perm)

    # --- Order tree branches ---
    for layer in k:-1:1 ###or 1?
        for tup in permutations(1:size(G,1), layer-1)
            if isa(tup, Int) #convert tuples of length 1 to lists
                tup = [tup]
            end
            tup2 = re_index(copy(tup))
            tstr = ""
            for i in tup2
                tstr = tstr * "[2]" * string([i]) * "[1]"
            end
            # tstr = "sort!(" * tstr * "[2])"
            # exec(tstr)
            tstr = tstr * "[2]"
            sort_in_place!(tree,tstr)
        end
    end
    return (key[2], tree)
end

# println(key_cospectrum([0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0],2))
Kriegman_1 = [0 0 0 1 0 1 0 1 0 1;
             0 0 0 0 1 1 0 1 0 1;
             0 0 0 0 0 1 1 1 1 0;
             1 0 0 0 0 0 1 0 1 1;
             0 1 0 0 0 0 1 0 1 1;
             1 1 1 0 0 0 0 1 0 0;
             0 0 1 1 1 0 0 0 1 0;
             1 1 1 0 0 1 0 0 0 0;
             0 0 1 1 1 0 1 0 0 0;
             1 1 0 1 1 0 0 0 0 0]

Kriegman_2 = [0 0 0 1 0 1 0 1 1 0;
             0 0 0 0 1 1 0 0 1 1;
             0 0 0 0 0 1 1 1 0 1;
             1 0 0 0 0 0 1 1 1 0;
             0 1 0 0 0 0 1 0 1 1;
             1 1 1 0 0 0 0 1 0 0;
             0 0 1 1 1 0 0 0 0 1;
             1 0 1 1 0 1 0 0 0 0;
             1 1 0 1 1 0 0 0 0 0;
             0 1 1 0 1 0 1 0 0 0]


# println(key_cospectrum(Kriegman_1,1))
# println(key_cospectrum(Kriegman_2,2))
# println(key_cospectrum(Kriegman_1,2)==key_cospectrum(Kriegman_2,2))