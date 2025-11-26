using Combinatorics
using LinearAlgebra

include("graph6_translation.jl")

"""
    inverse_permutation(p::Vector{Int}) -> Vector{Int}

Compute the inverse of a permutation p.
If p maps i -> p[i], then the inverse maps p[i] -> i.

Example: inverse_permutation([3, 1, 2]) = [2, 3, 1]
  because p[1]=3, p[2]=1, p[3]=2, so inverse[3]=1, inverse[1]=2, inverse[2]=3
"""
function inverse_permutation(p::Vector{Int})
    n = length(p)
    p_inv = zeros(Int, n)
    
    # For each position i, set p_inv[p[i]] = i
    # This creates the inverse mapping
    for i in 1:n
        p_inv[p[i]] = i
    end
    
    return p_inv
end

"""
    exec(string)

Legacy function that evaluates a string as Julia code.
Note: This uses eval() which can be unsafe. The codebase has moved away from
using this in favor of set_in_place! and get_in_place for safer nested access.
"""
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

"""
    spec_tuple(G::AbstractMatrix) -> Vector{Int}

Compute a signature tuple for graph G based on traces of matrix powers.
For an n×n adjacency matrix G, computes [tr(G²), tr(G³), ..., tr(Gⁿ)].

This is used as a graph invariant - graphs with different spectra will have
different spec_tuple values. Used in the cospectrality testing.
"""
function spec_tuple(G::AbstractMatrix)
    A = copy(G)
    out = Int[]
    # Compute traces of G², G³, ..., Gⁿ (where n = size(G,1))
    for i in 1:size(G, 1)-1
        A *= G  # A becomes G^(i+1)
        push!(out, tr(A))  # Store trace of G^(i+1)
    end
    return out
end

"""
    loop_diff(G, H) -> Int

Compute the squared difference between spec_tuples of two graphs G and H.
Returns sum of (spec_tuple(G)[i] - spec_tuple(H)[i])² for all i.

Used to measure how different two graphs are based on their spectra.
Smaller values indicate more similar graphs.
"""
function loop_diff(G, H)
    total = 0
    G_tup = spec_tuple(G)
    H_tup = spec_tuple(H)
    for i in eachindex(G_tup)
        total += (G_tup[i]-H_tup[i])^2
    end
    return total
end

"""
    empty_tree(n, depth) -> Vector

Create a nested tree structure used to store cospectrality information.
The tree structure is: [value, children] where:
  - value: an integer (initially n, later replaced with color indices)
  - children: a list of subtrees, each wrapped in a list [subtree]

Structure: Each node is [n, [[subtree1], [subtree2], ..., [subtree_n]]]
where each subtree has the same structure recursively.

This creates a tree where:
  - Root has value n and n children
  - Each child is a tree with value (n-1) and (n-1) children
  - Continues until depth reaches 0 or n reaches 0

Used in key_cospectrum() to store the canonical form of the graph's
cospectrality pattern.
"""
function empty_tree(n, depth)
    if n < 0 || depth < 0 || !(n == floor(n)) || !(depth == floor(depth))
        throw(ErrorException("Nonnegative integers only"))
    end
    # Base case: return just the value (leaf node)
    if n == 0 || depth == 0
        return n
    else
        # Recursive case: create n children, each containing a subtree
        children = []
        for j in 1:n
            # Each child is wrapped in a list: [subtree]
            push!(children, [empty_tree(n - 1, depth - 1)])
        end
        # Return [value, children] structure
        out = [n, children]
        return out
    end
end

"""
    re_index(lst::Vector{Int}) -> Vector{Int}

Re-index a permutation to allow for "renaming" that happens after vertex deletions

Used in key_cospectrum() to normalize vertex combinations before storing
them in the tree structure.
"""
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
        # This creates a canonical form where indices are "compressed"
        lst[i] -= count
    end
    return lst
end

"""
    relabel_tree(tree, vertices::Int, k::Int, bijection) -> tree

Apply a bijection (relabeling function) to all node values in the tree.
This is used during canonicalization to remap color indices according to
the sorted order of spectra.

The function:
1. Iterates through all layers (1 to k) and all vertex combinations
2. For each combination, finds the old color value
3. Applies the bijection to get the new color
4. Updates all permutations of that combination with the new color

The bijection maps old color indices to new ones based on the canonical
ordering of spectra (sorted by length and lexicographically).
"""
function relabel_tree(tree, vertices::Int, k::Int, bijection)
    # Process each layer of the tree
    for layer in 1:k
        # For each combination of 'layer' vertices
        for tup in combinations(1:vertices, layer)  # same as range(vertices) in Python
            # --- Find old colour ---
            # Normalize the combination using re_index
            altered_tup = re_index(copy(tup))
            # Build path string to navigate tree: "[2][p1][1][2][p2][1]..."
            # [2] = go to children, [p] = p-th child, [1] = get value
            tstr = ""
            for p in altered_tup
                tstr = tstr * "[2]" * string([p]) * "[1]"
            end
            # If not at deepest layer, add final [1] to access value
            if layer < k
                tstr = tstr * "[1]"
            end
            # Get the old color value at this location
            old_colour = get_in_place(tree,tstr)
            # Apply bijection to get new color
            new_colour = bijection[old_colour]

            #--- Relabel all permutations of tup ---
            # Update all permutations with the new color (to maintain symmetry)
            for perm in permutations(tup)
                altered_perm = re_index(copy(perm))
                tstr = ""
                for p in altered_perm
                    tstr = tstr * "[2]" * string([p]) * "[1]"
                end
                if layer < k
                    tstr = tstr * "[1]"
                end
                # Set the new color value at this location
                set_in_place!(tree, tstr, new_colour)
            end
        end
    end
    return tree
end

"""
    key_cospectrum(G, k::Int, test=spec_tuple) -> (Vector, tree)

Compute a "key" representation of graph G's cospectrum
up to depth k. This key is invariant under graph isomorphism and can be used
to test if two graphs are cospectral (have the same key).

The algorithm:
1. Builds a tree structure storing spectra of all vertex-deleted subgraphs
2. Assigns "colors" (indices) to each unique spectrum (stored in the key)
3. Canonically relabels colors based on sorted spectrum order
4. Orders tree branches to create a canonical form

Returns: (key[2], tree) where:
  - key[2]: sorted list of unique spectra (the "key")
  - tree: canonical tree structure encoding the cospectrality pattern

Two graphs are cospectral if and only if their keys are equal.

Parameters:
  - G: adjacency matrix of the graph
  - k: depth of cospectrality to compute (number of vertices to delete)
  - test: function to compute graph signature (default: spec_tuple)
"""
function key_cospectrum(G, k::Int, test=spec_tuple)
    # Base case: k=0 means no vertex deletion, just return the graph's spectrum
    if k == 0
        spec = test(G)
        return ([[spec], [1]])
    end

    # Initialize: create empty tree structure and key storage
    # tree structure: [value, children] where value will be color indices
    tree = empty_tree(size(G, 1), k)
    # key[1]: list of [spectrum, index] pairs for lookup
    # key[2]: list of unique spectra
    key = [[], []]

    # Process the full graph (no vertices deleted)
    stup = test(G)
    push!(key[1], [stup, 1])  # Store spectrum with index 1
    push!(key[2], stup)        # Add to unique spectra list
    tree[1] = 1                # Set root value to color 1

    # Process all layers: delete 1 vertex, then 2, ..., up to k vertices
    for layer in 1:k
        # For each combination of 'layer' vertices to delete
        for tup in combinations(1:size(G,1), layer)
            # Create submatrix by deleting rows/columns corresponding to tup
            M = G[1:end .∉ [tup], 1:end .∉ [tup]]
            stup = test(M)  # Compute spectrum of subgraph
            
            # If this is a new spectrum, add it to the key
            if !(stup in key[2])
                push!(key[1], [stup, length(key[2])+1])  # Store with next index
                push!(key[2], stup)                      # Add to unique list
            end

            # Find the color (index) assigned to this spectrum
            colour = findfirst(x -> x == stup, key[2])

            # Store this color for all permutations of the vertex combination
            # (to maintain symmetry - the order of deleted vertices shouldn't matter)
            for perm in permutations(tup)
                if isa(perm, Int) #convert tuples of length 1 to lists
                    perm = [perm]
                end
                # Build path string to navigate tree structure
                tstr = ""
                altered_perm = re_index(copy(perm))  # Normalize indices
                for p in altered_perm
                    # Path format: "[2][p][1]" means: go to children[2], then child p, then value[1]
                    tstr = tstr * "[2]" * string([p]) * "[1]"
                end
                # If not at deepest layer, add final [1] to access the value
                if layer < k #index???
                    tstr = tstr * "[1]"
                end
                # Store the color at this location in the tree
                set_in_place!(tree, tstr, colour)
            end
        end
    end

    # --- Canonical relabeling ---
    # Sort spectra to create a canonical ordering:
    # First by negative length (longer spectra first), then lexicographically
    tmp = deepcopy(key[2])
    sort!(tmp; by = x -> (-length(x),x))  #sort lex then by reverse length
    key[2] = tmp

    # Build bijection: map old color indices to new ones based on sorted order
    bijection = Int[]
    for i in 1:length(key[2]) #2:...???? instead
        # Find position of sorted spectrum[i] in original key[1]
        pos = findfirst(x -> x[1] == key[2][i], key[1])#???
        # Get the old index and add to bijection
        push!(bijection, key[1][pos][2])
    end

    # Apply inverse permutation to get mapping from old -> new colors
    bijection_perm = inverse_permutation(bijection)
    # Relabel entire tree with new canonical colors
    tree = relabel_tree(tree, size(G,1), k, bijection_perm)

    # --- Order tree branches ---
    # Sort children at each node to create fully canonical form
    # Process from deepest layer to shallowest (k down to 1)
    for layer in k:-1:1 ###or 1?
        # For each combination of (layer-1) vertices
        for tup in permutations(1:size(G,1), layer-1)
            if isa(tup, Int) #convert tuples of length 1 to lists
                tup = [tup]
            end
            # Build path to the node whose children we want to sort
            tup2 = re_index(copy(tup))
            tstr = ""
            for i in tup2
                tstr = tstr * "[2]" * string([i]) * "[1]"
            end
            # Append [2] to path to access children list for sorting
            tstr = tstr * "[2]"
            # Sort children by their values (colors) to create canonical order
            sort_in_place!(tree,tstr)
        end
    end
    return (key[2], tree)
end

# ============================================================================
# TEST CASES
# ============================================================================

# Example: Test on a simple 4-vertex graph
# println(key_cospectrum([0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0],2))

# Kriegman graphs: Two 10-vertex graphs that are 1-cospectral, found by Aaron Kriegman (2024)
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

# Test cases (commented out):
# println(key_cospectrum(Kriegman_1,1))
# println(key_cospectrum(Kriegman_2,2))
# Check if graphs are cospectral (should return true for cospectral graphs):
# println(key_cospectrum(Kriegman_1,2)==key_cospectrum(Kriegman_2,2))