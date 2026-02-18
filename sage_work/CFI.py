from sage.all import *
from cospectrality import *
from copy import copy
from key_cospectrum import *
import time

def get_block_bounds(blocks, i, j):
    row_start = sum(block.nrows() for block in blocks[:i])
    row_end = row_start + blocks[i].nrows()
    
    col_start = sum(block.ncols() for block in blocks[:j])
    col_end = col_start + blocks[j].ncols()
    
    return row_start, row_end, col_start, col_end

def set_block(M, blocks, i, j, new_block):
    #Set block (i, j) in M to new_block (must match shape).
    r0, r1, c0, c1 = get_block_bounds(blocks, i, j)
    if new_block.nrows() != r1 - r0 or new_block.ncols() != c1 - c0:
        raise ValueError(f"Shape mismatch: expected ({r1 - r0}, {c1 - c0})")
    M[r0:r1, c0:c1] = new_block

def gadget(d):
    M_d0 = [bin(i)[2:] for i in range(2**(d-1))]
    M_d = []
    for i in M_d0:
        s = sum([int(digit) for digit in i])%2
        string = i + str(s)
        while len(string) < d:
            string = '0' + string
        M_d.append(string)

    F_d = Matrix([[0 for i in range(2*d+len(M_d))] for j in range(2*d+len(M_d))])
    for i in range(d):
        for j in range(2**(d-1)):
            if M_d[j][i] == '0':
                F_d[2*i,2*d+j] = 1
            else:
                F_d[2*i+1,2*d+j] = 1
    
    F_d = F_d + F_d.transpose()

    return(F_d)

def CFI(G, twist = 0):
    #find twist point if we will twist
    i = 0
    j = 0
    twist_point = [i,j]
    while G[i,j] == 0:
        if j == G.nrows() - 1:
            j = 0
            i = i + 1
        else:
            j = j + 1
        twist_point = [i,j]

    #now make all necessary gadgets    
    degs = sorted(list(set(Graph(G).degree_sequence())))
    gadgets = []
    for i in range(degs[-1]+1):
        if i in degs:
            gadgets.append(gadget(i))
        else:
            gadgets.append(Matrix([[0]]))
    
    #compute their sizes for when we make block matrix
    gadget_sizes = [m.nrows() for m in gadgets]
    degs = [sum(G[i]) for i in range(G.nrows())]
    total_size = 0
    for d in degs:
        total_size = total_size + gadget_sizes[d]    

    blocks = [gadgets[degs[i]] for i in range(G.nrows())]

    M = block_diagonal_matrix(blocks)
    #set the intersection stuff to join a's and b's
    for i in range(G.nrows()):
        for j in range( G.nrows()):
                if G[i,j] == 1  :
                    r0, r1, c0, c1 = get_block_bounds(blocks, i, j)
                    dims = [r1-r0,c1-c0]
                    new_block = Matrix([[0 for j in range(dims[1])] for j in range(dims[0])])
                    #find out which a_i joins which a_i
                    row_a = sum(G[i][0:j])
                    col_a = sum(G[j][0:i])
                    if twist == 1 and ((i == twist_point[0] and j == twist_point[1]) or (j == twist_point[0] and i == twist_point[1])):
                        new_block[2*row_a,2*col_a+1] = 1
                        new_block[2*row_a+1,2*col_a] = 1
                    else:
                        new_block[2*row_a,2*col_a] = 1
                        new_block[2*row_a+1,2*col_a+1] = 1

                    set_block(M, blocks, i, j, new_block)
     
    return(M)

S4 = graphs.CompleteBipartiteGraph(1,4).adjacency_matrix()
S5 = graphs.CompleteBipartiteGraph(1,5).adjacency_matrix()
S6 = graphs.CompleteBipartiteGraph(1,6).adjacency_matrix()

K = [0,0] + [graphs.CompleteGraph(i).adjacency_matrix() for i in range(2,9)]


def complete_cfi_1_less(m, k):
    a = CFI(K[m],0)
    b = CFI(K[m],1)
    n = a.nrows()

    for i in range(2*(n/m) + 1):
        a_1 = a.delete_rows([i]).delete_columns([i])
        b_1 = b.delete_rows([i]).delete_columns([i])
        if key_cospectrum(a_1,k-1) != key_cospectrum(b_1,k-1):
            return(i,"is the problem")
        else:
            print(str(i)+"/"+str(n-1))
    print(str(k) + "-cospectral")
    return(True)

def list_del_tester(g1,g2,lst1,lst2,k):
    a_1 = g1.delete_rows(lst1).delete_columns(lst1)
    b_1 = g2.delete_rows(lst2).delete_columns(lst2)
    return(key_cospectrum(a_1,k) == key_cospectrum(b_1,k))


