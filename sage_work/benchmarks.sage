import os
from cospectrality import *
from sage.all import *
from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend

def get_filenames_from_folder(folder_path):
    # Get all filenames in the specified folder
    filenames = os.listdir(folder_path)
    
    # Return the list of filenames
    return filenames

def dimacs_to_adjacency_matrix(input):
    # Step 1: Parse the number of vertices (from the 'p' line)
    setup = input[0].split()
    vertices = int(setup[2])
    
    # Step 2: Create an empty graph with 'vertices' vertices
    G = Graph({0: [i for i in range(1,vertices+1)]} )
    
    # Step 3: Add edges to the graph based on the 'e' lines
    for edge in input:
        if edge.startswith('e'):
            # Extract vertex pairs from 'e' lines, e.g., 'e 1 3' -> (1, 3)
            parts = edge.split()
            vertex1 = int(parts[1])
            vertex2 = int(parts[2])
            G.add_edge(vertex1, vertex2)
    
    # Step 4: Get the adjacency matrix of the graph
    adj_matrix = G.adjacency_matrix()
    
    # Return the adjacency matrix
    return adj_matrix


wd = os.getcwd()
path = wd+"/cfi-rigid-z2"
files = (get_filenames_from_folder(path))

# for index in range(int(len(files)/2)):
for index in range(1):

    g1 = dimacs_to_adjacency_matrix(open(path+"/"+files[2*index], 'r').readlines())
    g2 = dimacs_to_adjacency_matrix(open(path+"/"+files[2*index + 1], 'r').readlines())

    print("graphs "+str(2*index)+" and "+str(2*index +1))
    print(str(g1.nrows())+ " vertices")
    print("lalala",multiplicity_k_cospectrum(g1,3)==multiplicity_k_cospectrum(g2,3))
    test = Graph(g1).is_isomorphic(Graph(g2))
    if test == True:
        print("isomorphic")
    # elif multiplicity_k_cospectrum(g1,2)!=multiplicity_k_cospectrum(g2,2):
    #     print("Not 2-cospectral")
    elif multiplicity_k_cospectrum(g1,3)!=multiplicity_k_cospectrum(g2,3):
        print("Not 3-cospectral")
    elif multiplicity_k_cospectrum(g1,4)!=multiplicity_k_cospectrum(g2,4):
        print("Not 4-cospectral")
    else:
        print(how_cospectral(g1,g2,1))