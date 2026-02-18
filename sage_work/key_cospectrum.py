from sage.all import *
from cospectrality import *
from copy import *

def empty_tree(n, depth):
    if n < 0 or depth < 0 or n != int(n) or depth != int(depth):
        raise TypeError("Nonnegative integers only")
    if n == 0 or depth == 0:
        return(n)
    else:
        children = []
        for j in range(n):
            children.append([empty_tree(n - 1, depth - 1)])
        out = [n, children]################
    return(out)

def relabel_tree(tree, vertices, k, bijection):
    for layer in range(1,k+1):
        for tup in Combinations([n for n in range(vertices)], Integer(layer)):
            #find old_colour
            altered_tup = list(copy(tup))
            altered_tup = re_index(altered_tup)
            tstr = "tree"
            for p in altered_tup:
                tstr = tstr + "[1]" + str([p]) + "[0]" # +1 due to index matching up
            if layer < k:
                tstr = tstr + "[0]"
            old_colour = eval(tstr)
            new_colour = bijection(old_colour)

            for perm in Permutations(tup):
                tstr = "tree"
                # account for shifting of indices after vertex deletions
                altered_perm = list(copy(perm))
                altered_perm = re_index(altered_perm)         

                for p in altered_perm:
                    tstr = tstr + "[1]" + str([p]) + "[0]" # +1 due to index matching up
                if layer < k:
                    tstr = tstr + "[0]"
                exec(tstr + "= new_colour")
    return(tree)

def key_cospectrum(G, k):
    if k == 0:
        return([[[spec_tuple(G),0]],[spec_tuple(G)]],[0])
    tree = empty_tree(G.nrows(), k)
    key = [[],[]]

    stup = spec_tuple(G)
    key[0].append([stup,0])
    key[1].append(stup)
    tree[0] = 0

    for layer in range(1,k+1):
        for tup in Combinations([n for n in range(G.nrows())], Integer(layer)):
            M = G.delete_rows(tup).delete_columns(tup)
            stup = spec_tuple(M)
            if stup not in key[1]:
                key[0].append([stup,len(key[1])])
                key[1].append(stup)
            position = 0
            while stup != key[1][position]:
                position = position + 1
            # colour = key[0][position][1]
            colour = position

            for perm in Permutations(tup):
                tstr = "tree"
                # account for shifting of indices after vertex deletions
                altered_perm = list(copy(perm))
                altered_perm = re_index(altered_perm)         

                for p in altered_perm:
                    tstr = tstr + "[1]" + str([p]) + "[0]" # +1 due to index matching up
                if layer < k:
                    tstr = tstr + "[0]"
                exec(tstr + "= colour")

    #now relabel the tree canonically
    tmp = key[1]
    tmp.sort() # lex order
    tmp.sort(key=len, reverse=True) #then reverse length order
    key[1] = tmp

    bijection = []
    for i in range(1,len(key[1])):
        pos = 0
        while key[0][pos][0] != key[1][i]:
            pos = pos + 1
        bijection.append(key[0][pos][1])
    bijection = Permutation(bijection).inverse()
    tree = relabel_tree(tree, G.nrows(), k, bijection)

    #now order tree branches
    for layer in range(k, 0, -1):
        for tup in Permutations([n for n in range(G.nrows())],Integer(layer-1)):#####
            tup2 = copy(tup)
            tup2 = re_index(tup2)
            tstr = "tree"
            for i in tup2:
                tstr = tstr + "[1]" + str([i]) + "[0]"
            tstr = tstr + "[1]"
            exec(tstr + ".sort()")
    

    key = key[1]
    return(key, tree)

def key_how_cospectral(G, H):
    switch = True
    i = 0
    while switch == True:
        switch = (key_cospectrum(G,i) == key_cospectrum(H,i))
        i = i+1
    return(i-2)