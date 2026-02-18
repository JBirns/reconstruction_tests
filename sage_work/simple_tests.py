from sage.all import *
from cospectrality import *
from copy import copy
from key_cospectrum import *

def edge_count(G):
    return(sum(G.list()))

def empty(G):
    if G.is_zero() == True:
        return(1)
    else:
        return(0)
    
def degseq(G):
    degs = []
    for i in range(G.nrows()):
        degs.append(sum(G[i]))
    tmp = sorted(degs)
    tmp.reverse()
    return(tmp)

def ds_1(G):
    n = G.nrows()
    output = []
    seqs = [] 

    for i in range(n):
        newseq = []
        for j in range(n):
            if G[i,j] == 1:
                newseq.append(sum(G[j]))
        newseq = sorted(newseq)
        seqs.append(newseq)

    seqs = sorted(seqs)
    output.append([seqs[0],1])
    for i in range(1,len(seqs)):
        if seqs[i]==seqs[i-1]:
            output[-1][1] = output[-1][1] + 1
        else:
            output.append([seqs[i],1])

    # output = sorted(output)
    return(output)

def edge_3(G):
    n = G.nrows()
    output = []
    ds1s = [ds_1(G.delete_rows([i]).delete_columns([i])) for i in range(n)]
    ds1s = sorted(ds1s)
    output.append([ds1s[0],1])
    for i in range(1,len(ds1s)):
        if ds1s[i]==ds1s[i-1]:
            output[-1][1] = output[-1][1] + 1
        else:
            output.append([ds1s[i],1])

    # output = sorted(output)
    return(output)

def edge_4(G):
    n = G.nrows()
    output = []
    ds1s = [edge_3(G.delete_rows([i]).delete_columns([i])) for i in range(n)]
    ds1s = sorted(ds1s)
    output.append([ds1s[0],1])
    for i in range(1,len(ds1s)):
        if ds1s[i]==ds1s[i-1]:
            output[-1][1] = output[-1][1] + 1
        else:
            output.append([ds1s[i],1])

    # output = sorted(output)
    return(output)

def simple_tree(G, k, test):
    tree_base = []
    n = G.nrows()

    if k == 0:
        return([test(G)])
    tree = empty_tree(G.nrows(), k)

    for tup in Combinations([t for t in range(n)], Integer(k)):
        M = G.delete_rows(tup).delete_columns(tup)
        stup = test(M)

        for perm in Permutations(tup):
            tstr = "tree"
            # account for shifting of indices after vertex deletions
            altered_perm = list(copy(perm))
            altered_perm = re_index(altered_perm)         

            for p in altered_perm:
                tstr = tstr + "[1]" + str([p]) + "[0]" # +1 due to index matching up
            exec(tstr + "= stup")

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

    return(tree)

def simple_test(G,H,k,test):
    if simple_tree(G,k,test) == simple_tree(H,k,test):
        return(True)
    else:
        return(False)