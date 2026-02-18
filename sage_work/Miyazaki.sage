from sage.all import *
from cospectrality import *

Miyazaki_I = [[0, 4], [1, 4], [2, 5], [3, 5], [0, 6], [2, 6], [1, 7], [3, 7], [0, 8], [3, 8], [1, 9], [2, 9], [4, 6], [5, 7], [10, 14], [11, 14], [12, 15], [13, 15], [10, 16], [12, 16], [11, 17], [13, 17], [10, 18], [13, 18], [11, 19], [12, 19], [14, 16], [15, 17], [8, 18], [9, 19]]
Miyazaki_II = [[0, 4], [1, 4], [2, 5], [3, 5], [0, 6], [2, 6], [1, 7], [3, 7], [0, 8], [3, 8], [1, 9], [2, 9], [4, 6], [5, 7], [10, 14], [11, 14], [12, 15], [13, 15], [10, 16], [12, 16], [11, 17], [13, 17], [10, 18], [13, 18], [11, 19], [12, 19], [14, 16], [15, 17], [8, 19], [9, 18]]

def czm_to_sage_mat(edges):

    num_vertices = 0

    for e in edges:
        if e[1]+1 > num_vertices:
            num_vertices = e[1]+1

    adj_matrix = Matrix(ZZ, num_vertices, num_vertices, 0)
    
    # For each edge in the list, set the corresponding entries to 1 (assuming an undirected graph)
    for edge in edges:
        u, v = edge
        adj_matrix[u, v] = 1
        adj_matrix[v, u] = 1  # Since the graph is undirected
    
    return adj_matrix

MI = czm_to_sage_mat(Miyazaki_I)
MII = czm_to_sage_mat(Miyazaki_II)

print(how_cospectral(MI,MII))