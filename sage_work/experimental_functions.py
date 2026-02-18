from sage.all import *
from cospectrality import *
from copy import copy
from key_cospectrum import *
from tree_compressor import *
from sage.graphs.trees import TreeIterator

def structural_comparison(Mats, k):
    for m in range(len(Mats)):
        a = tree_structure(Mats[m],k)
        for n in range(m+1, len(Mats)):
            b = tree_structure(Mats[n],k)
            if b == a:
                print("There exist a pair with the same tree structure")
            print(m+1, n+1, len(Mats))

def subgraphs_count(G,k):
    n = G.nrows()
    specs = []
    counts = []
    for tup in Combinations([m for m in range(n)], Integer(n-k)):
        M = G.delete_rows(tup).delete_columns(tup)
        stup = spec_tuple(M)
        if stup not in specs:
            specs.append(stup)
            counts.append(1)
        else:
            i = 0
            while specs[i] != stup:
                i = i+1
            counts[i] == counts[i] + 1
    return(specs,counts)

def subgraphs_compare(Mats, k):
    for i in range(len(Mats)):
        a = subgraphs_count(Mats[i], k)
        for j in range(i+1,len(Mats)):  
            b = subgraphs_count(Mats[j],k)
            if a == b:
                print("Match!")
            else:
                print(i+1,j+1,len(Mats))

def question_2_3_checker(G):
    n = G.nrows()
    combs = Combinations([m for m in range(n)], 4)
    specs = []
    coloured_subgs = []
    for tup in combs:
        M = G.delete_rows(tup).delete_columns(tup)
        spec = spec_tuple(M)
        if spec not in specs:
            specs.append(spec)
        coloured_subgs.append(specs.index(spec))
    
    for i in range(len(combs)):
        print(i,len(combs))
        M1 = G.delete_rows(combs[i]).delete_columns(combs[i])
        M1 = Graph(M1)
        tup_i_2 = [combs[i][0], combs[i][1]]
        tup_i_3 = [combs[i][0], combs[i][1], combs[i][2]]
        mitest2 = G.delete_rows(tup_i_2).delete_columns(tup_i_2)
        mitest3 = G.delete_rows(tup_i_3).delete_columns(tup_i_3)
        ki2 = key_cospectrum(mitest2,2)
        ki3 = key_cospectrum(mitest3,1)
        for j in range(i+1,len(combs)):
            tup_j_2 = [combs[j][0], combs[j][1]]
            tup_j_3 = [combs[j][0], combs[j][1], combs[j][2]]
            mjtest2 = G.delete_rows(tup_j_2).delete_columns(tup_j_2)
            mjtest3 = G.delete_rows(tup_j_3).delete_columns(tup_j_3)
            if coloured_subgs[i] == coloured_subgs[j] and key_cospectrum(mjtest2,2) == ki2 and key_cospectrum(mjtest3,1) == ki3:
                M2 = G.delete_rows(combs[j]).delete_columns(combs[j])
                M2 = Graph(M2)
                if not M1.is_isomorphic(M2):
                    return(False)
    return(True)

def check_trees(trees):
    max = -1
    for i in range(len(trees)):
        for j in range(i+1, len(trees)):
            val = key_how_cospectral(trees[i],trees[j])
            if val > max:
                max = val
                print(val,i,j)
    return(max)

def children_level_k(G, k):
    if k == 0:
        return(spec_tuple(G))
    
    children = []
    
    for tup in Combinations([n for n in range(G.nrows())], Integer(k)):
            M = G.delete_rows(tup).delete_columns(tup)
            stup = spec_tuple(M)
            if stup not in children:
                children.append(stup)

    children.sort()
    
    return(children)
    
def children_intersection_checker(Mats, k):
    for i in range(len(Mats)):
        A = children_level_k(Mats[i], k)
        for j in range(i+1, len(Mats)):
            B = children_level_k(Mats[j], k)
            for s in B:
                if s in A:
                    return(False)
            print(i+1,j+1,len(Mats))
    return(True)


def coalesce_trees(A,B):
    # Get dimensions of the input matrices
    n = A.nrows()
    m = B.nrows()
    # Ensure both matrices are square
    if A.ncols() != n or B.ncols() != m:
        raise ValueError("Both matrices must be square.")
    result = Matrix(ZZ, n + m - 1, n + m - 1, 0)
    result[0:n, 0:n] = A[0:n, 0:n]
    result[n-1:n-1+m, n-1:n-1+m] = B    
    return result

def RS_trees(trees):
    s_trees = []
    r_trees = []
    R = Matrix([[0,1,0,0,0,0,0,0,0],
                [1,0,1,0,0,0,0,0,0],
                [0,1,0,0,0,0,0,0,1],
                [0,0,0,0,1,0,0,0,1],
                [0,0,0,1,0,1,0,1,0],
                [0,0,0,0,1,0,1,0,0],
                [0,0,0,0,0,1,0,0,0],
                [0,0,0,0,1,0,0,0,0],
                [0,0,1,1,0,0,0,0,0]])
    
    S = Matrix([[0,1,0,0,0,0,0,0,0],
                [1,0,1,0,0,0,0,0,0],
                [0,1,0,1,0,0,0,0,0],
                [0,0,1,0,1,0,0,0,0],
                [0,0,0,1,0,1,0,0,0],
                [0,0,0,0,1,0,0,1,1],
                [0,0,0,0,0,0,0,0,1],
                [0,0,0,0,0,1,0,0,0],
                [0,0,0,0,0,1,1,0,0]])
    for t in trees:
        #attach s or r to vertex 1 of t
        s_trees.append(coalesce_trees(S,t))
        r_trees.append(coalesce_trees(R,t))
    return(r_trees,s_trees)

def check_cospectral_trees(n, brakes):
    max = -1
    trees = []
    for t in TreeIterator(n):
        t = t.adjacency_matrix()
        trees.append(t)
        if len(trees) == brakes:
            break
    pairs = RS_trees(trees)
    for i in range(len(pairs[0])-1):
        val = key_how_cospectral(pairs[0][i],pairs[1][i])        
        if val > max:
            max = val
        print(val , max , i+1)
    return(max)




