from sage.all import *
from copy import copy

def spec_tuple(G): #G an adj. matrix
    A = copy(G)
    out = []
    for i in range(1,G.nrows()):
        A = A * G
        out.append(A.trace())
    return(out)        


def childbear(G):
    return([spec_tuple(G),
            [spec_tuple(G.delete_rows([i]).delete_columns([i])) for i in range(G.nrows())]])


def re_index(lst):
    # Traverse the list from right to left
    for i in range(len(lst) - 1, -1, -1):
        # Count how many elements to the right are less than lst[i]
        count = sum(1 for j in range(i) if lst[j] < lst[i])
        # Decrease the current term by the count
        lst[i] = lst[i] - count
    return(lst)

def generationbear(G, old_cospectrum):
    k = 0
    cstr = "old_cospectrum" + "[0]"
    while not (isinstance(eval(cstr), Integer) or isinstance(eval(cstr), int)):
        cstr = cstr[:-3]
        cstr = cstr + "[1][0][0]"
        k = k + 1

    if k>0:
        for tup in Permutations([n for n in range(G.nrows())],Integer(k)):
            tup2 = copy(tup)
            tup2 = re_index(tup2)
            cstr = "old_cospectrum"
            for i in range(len(tup2)):
                cstr = cstr + "[1]" + str([tup2[i]])
            exec(cstr + "=" + str(childbear(G.delete_rows(tup).delete_columns(tup))))    
    else:
        old_cospectrum = childbear(G)
    
    return(old_cospectrum)


def order_cospec(G, cospectrum, k):
    for layer in range(k, 0, -1): #layer was i and layer2 was k
        for tup in Permutations([n for n in range(G.nrows())],Integer(layer-1)):#####
            tup2 = copy(tup)
            tup2 = re_index(tup2)
            cstr = "cospectrum[1]"
            for i in range(len(tup2)):
                cstr = cstr + str([tup2[i]]) + "[1]"
            exec(cstr + ".sort()")
    return(cospectrum)


def cospectrum(G, k):
    if k > G.nrows()-4:
        raise TypeError("k excessively large")
    cospectrum = spec_tuple(G)
    if k == 0:
        return(cospectrum)
    else:
        for i in range(k):
            cospectrum = generationbear(G, cospectrum)

    #now order branches
    for i in range(k,0,-1):
        order_cospec(G,cospectrum,i)
        
    return(cospectrum)


def cospectrality_test(G, H, k):
    return(cospectrum(G, k) == cospectrum(H, k))

def how_cospectral(G, H, counter=0):
    k = 0
    co_G = spec_tuple(G)
    co_H = spec_tuple(H)
    if co_G != co_H:
        return("Not cospectral")
    else:
        while cospectrum(G,k) == cospectrum(H,k):
            if counter == 1:
                print(str(k)+"-cospectral")
            k = k + 1
            # co_G = generationbear(G, co_G)
            # co_H = generationbear(H, co_H)
            # for i in range(k,0,-1):
            #     order_cospec(G,co_G,i)
            #     order_cospec(H,co_H,i)
            if k == G.nrows() - 3:
                return("Isomorphic")
        return(k-1)


def multiplicity_k_cospectrum(G, k):

    counts = []
    spec_tuples = []

    for tup in Combinations([n for n in range(G.nrows())],Integer(k)):
        M = G.delete_rows(tup).delete_columns(tup)
        # stup = spec_tuple(M)
        stup = M.eivenvalues()
        if stup in spec_tuples:
            i = 'temp'
            for j in range(len(spec_tuples)):
                if stup == spec_tuples[j]:
                    i = j
                    break
            tmp = copy(counts[i])
            counts[i] = tmp + 1
        else:
            spec_tuples.append(stup)
            counts.append(1)
        
    output = [[counts[i],spec_tuples[i]] for i in range(len(counts))]
    output.sort()
    return(output)