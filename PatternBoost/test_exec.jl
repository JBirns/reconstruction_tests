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

ar = [[2], [9,4,1], [9,3,5,5]]
# sort!(ar)
sort_in_place!(ar,"")
get_in_place(ar,"[2][2]")