


# Function to add the transpose of the matrix to itself
function addTranspose!(adjMatrix)
    for j in 1:length(adjMatrix)
        for i in 1:j-1
            adjMatrix[i][j] = adjMatrix[j][i]
        end
    end
end

function graph6ToAdj(graph6::AbstractString)
    # vertices + 63 = first char
    vertices = Int(Char(graph6[1])) - 63

    bin_list = ""

    # Turn into 6-bit pieces
    for i in 2:length(graph6)
        bin_list *= lpad(string(Int(Char(graph6[i])) - 63, base=2), 6, '0')
    end

    adjMatrix = []

    # Unpad on right until we have bottom-left diagonal
    num_in_bot_left_diag = 0
    for i in 1:vertices-1
        num_in_bot_left_diag += i
    end

    bot_left_diag = bin_list[1:num_in_bot_left_diag]

    idx = 1
    for i in 1:vertices
        sub_adjMatrix = zeros(Int, vertices)
        for j in 1:i-1
            sub_adjMatrix[j] = parse(Int, bin_list[idx])
            idx += 1
        end
        push!(adjMatrix, sub_adjMatrix)
    end

    addTranspose!(adjMatrix)

    return adjMatrix
end

function graph6ToAdj(graph6::AbstractString)
    # Calculate the number of vertices
    vertices = Int(Char(graph6[1])) - 63

    bin_list = ""

    # Turn the Graph6 string into a list of 6-bit binary pieces
    for i in 2:length(graph6)
        bin_list *= lpad(string(Int(Char(graph6[i])) - 63, base=2), 6, '0')
    end

    # Initialize the adjacency matrix as a zero matrix of size vertices x vertices
    adjMatrix = zeros(Abstract, vertices, vertices)

    # Unpad on the right until we have the lower-left diagonal
    num_in_bot_left_diag = 0
    for i in 1:vertices-1
        num_in_bot_left_diag += i
    end

    bot_left_diag = bin_list[1:num_in_bot_left_diag]

    idx = 1
    for i in 1:vertices
        for j in 1:i-1  # Only fill the lower triangular part
            adjMatrix[i, j] = parse(Int, bot_left_diag[idx])
            idx += 1
        end
    end

    # Fill the upper triangular part (transpose)
    for i in 1:vertices
        for j in 1:i-1
            adjMatrix[j, i] = adjMatrix[i, j]
        end
    end

    return vertices, adjMatrix
end