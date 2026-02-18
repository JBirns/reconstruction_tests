from cospectrality import *
from key_cospectrum import *
from sage.all import *
from random import randint
from CFI import *
import time

def random_cosp_test(G, H, k):
    if spec_tuple(G) != spec_tuple(H):
        return("Not cospectral")
    
    n = G.nrows()
    seq = []
    while len(seq) < k:
        tmp = randint(0,n-1)
        if tmp not in seq:
            seq.append(tmp)
    
    specs_G = [spec_tuple(G.delete_rows(seq[0:i]).delete_columns(seq[0:i])
                         ) for i in range(len(seq)+1)]
    
    legit_tups = []
    for i in range(1,k+1):
        for tup in Permutations([t for t in range(n)],Integer(i)):
            if i == 1 and spec_tuple(
                H.delete_rows(tup).delete_columns(tup)) in specs_G:
                legit_tups.append(list(tup))
            elif tup[0:i-1] in legit_tups and spec_tuple(
                H.delete_rows(tup).delete_columns(tup)) in specs_G:
                legit_tups.append(list(tup))
        if legit_tups == []:
            return(False, "Not 1-cospectral")
        elif len(legit_tups[-1]) < i:
            return(False, "Not "+str(i)+"-cospectral")
    
    return(True, "May be "+str(k)+"-cospectral")

# for i in range(10):
#     start = time.time()
#     print(random_cosp_test(
#         CFI(S5,0),CFI(S5,1),
#         4))
#     end = time.time()
#     print(end - start)
#     print("test",i)
