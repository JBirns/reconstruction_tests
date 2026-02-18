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

# print(key_how_cospectral(CFI(S4,0),CFI(S4,1)))
# print(key_how_cospectral(CFI(S5,0),CFI(S5,1)))
# print(key_how_cospectral(CFI(S6,0),CFI(S6,1)))

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
    
def is_weakly_decreasing(lst):
    return all(x >= y for x, y in zip(lst, lst[1:]))

def is_weakly_decreasing_nonzero(lst):
    nonzeros = [x for x in lst if x != 0]
    return all(x >= y for x, y in zip(nonzeros, nonzeros[1:]))


m,k = 5,0
a,b = CFI(K[m],0),CFI(K[m],1)
n = a.nrows()
switch = 0
Fd = n/m
Od = 2*(m-1)
Md = Fd-Od

# for i in range(n):
#     if not list_del_tester(a,b,[Fd+Od,i],[Fd+Od,i],k):
#         for j in range(n):
#             if list_del_tester(a,b,[Fd+Od,i],[Fd+Od,j],k):
#                 print(i,j)
#     else:
#         print(i,"is fine")

# for i in range(6,10):
#     for j in range(10,n-1):
#         if not list_del_tester(a,b,[i,j],[i,j+1],k-1):
#             print(i,j,"is the problem")
#             switch = 1
#             break
#         else:
#             print(str(i)+"/"+str(n-1))
#     if switch == 1:
#         break


def complete_CFI_explicit_bijection(m, dels):
    if dels < 1:
        raise TypeError("dels<1")
    a,b = CFI(K[m],0),CFI(K[m],1)
    n = a.nrows()
    Fd = n/m
    Od = 2*(m-1)
    Md = Fd-Od

    if dels == 1:
        return([[i] for i in range(n)],[[i] for i in range(n)])
    
    pair_list = []
    for t in Permutations([t for t in range(a.nrows())],Integer(2)):
        pair_list.append(list(t))
    pair_bijection = []
    for pair in pair_list:
        i,j = pair[0],pair[1]

        #i < j
        if i in range(0,2) and j == Fd:
            new_pair = [i,j+1]
        elif i in range(0,2) and j == Fd+1:
            new_pair = [i,j-1]
        elif i in range(0,2) and j in range(Fd+Od,Fd+Od+Md/2): #is this pairs or cycle (gone with cycle)?
            new_pair = [i,j+Md/2]
        elif i in range(0,2) and j in range(Fd+Od+Md/2,2*Fd):
            new_pair = [i,j-Md/2]

        elif i in range(Od,Fd) and j == Fd:
            new_pair = [i,j+1]
        elif i in range(Od,Fd) and j == Fd+1:
            new_pair = [i,j-1]
        elif i in range(Od,Fd) and j in range(Fd+Od,Fd+Od+Md/2): #is this pairs or cycle (gone with cycle)?
            new_pair = [i,j+Md/2]
        elif i in range(Od,Fd) and j in range(Fd+Od+Md/2,2*Fd):
            new_pair = [i,j-Md/2]

        #i > j
        elif i in range(Fd,Fd+2) and j in range(0,2):
            new_pair = [i, (j+1)%2]
        elif i in range(Fd+Od,2*Fd) and j in range(0,2):
            new_pair = [i, (j+1)%2]
        
        elif i in range(Fd, Fd+2) and j in range(Od,Od+Md/2):
            new_pair = [i,j+Md/2]
        elif i in range(Fd, Fd+2) and j in range(Od+Md/2,Fd):
            new_pair = [i,j-Md/2]
        elif i in range(Fd+Od,2*Fd) and j in range(Od,Od+Md/2):
            new_pair = [i,j+Md/2]
        elif i in range(Fd+Od,2*Fd) and j in range(Od+Md/2,Fd):
            new_pair = [i,j-Md/2]

        else:
            new_pair = [i,j]

        pair_bijection.append(new_pair)

    list_1 = []
    for t in Permutations([t for t in range(a.nrows())],Integer(dels)):
        list_1.append(list(t))
    list_2 = []

    for i in range(len(list_1)):
        tmp = copy(list_1[i])
        pair = tmp[0:2]
        del tmp[0:2]
        bij_pair = pair_bijection[pair_list.index(pair)]
        list_2.append(bij_pair + tmp)

    return(list_1,list_2)


def complete_CFI_explicit_bijection_2(m, dels):
    if dels < 1:
        raise TypeError("dels<1")
    a,b = CFI(K[m],0),CFI(K[m],1)
    n = a.nrows()
    Fd = n/m
    Od = 2*(m-1)
    Md = Fd-Od

    if dels == 1:
        return([[i] for i in range(n)],[[i] for i in range(n)])
    
    list_1 = []
    for t in Permutations([t for t in range(a.nrows())],Integer(dels)):
        list_1.append(list(t))
    list_2 = []

    bad_pairs = []
    good_pairs = []
    pair_list = []
    for t in Permutations([t for t in range(a.nrows())],Integer(2)):
        pair_list.append(list(t))
    
    for pair in pair_list:
        i,j = pair[0],pair[1]

        if i in range(0,2) and j == Fd:
            new_pair = [i,j+1]
        elif i in range(0,2) and j == Fd+1:
            new_pair = [i,j-1]
        elif i in range(0,2) and j in range(Fd+Od,Fd+Od+Md/2): #is this pairs or cycle (gone with cycle)?
            new_pair = [i,j+Md/2]
        elif i in range(0,2) and j in range(Fd+Od+Md/2,2*Fd):
            new_pair = [i,j-Md/2]

        elif i in range(Od,Fd) and j == Fd:
            new_pair = [i,j+1]
        elif i in range(Od,Fd) and j == Fd+1:
            new_pair = [i,j-1]
        elif i in range(Od,Fd) and j in range(Fd+Od,Fd+Od+Md/2): #is this pairs or cycle (gone with cycle)?
            new_pair = [i,j+Md/2]
        elif i in range(Od,Fd) and j in range(Fd+Od+Md/2,2*Fd):
            new_pair = [i,j-Md/2]

        #i > j
        elif i in range(Fd,Fd+2) and j in range(0,2):
            new_pair = [i, (j+1)%2]
        elif i in range(Fd+Od,2*Fd) and j in range(0,2):
            new_pair = [i, (j+1)%2]
        
        # elif i in range(Fd, Fd+2) and j in range(Od,Od+Md):####
        #     new_pair = [i,Fd+Od-1-j]
        elif i in range(Fd, Fd+2) and j in range(Od,Od+Md/2):#####
            new_pair = [i,j+Md/2]
        elif i in range(Fd, Fd+2) and j in range(Od+Md/2,Fd):#####
            new_pair = [i,j-Md/2]
        elif i in range(Fd+Od,2*Fd) and j in range(Od,Od+Md/2):
            new_pair = [i,j+Md/2]
        elif i in range(Fd+Od,2*Fd) and j in range(Od+Md/2,Fd):
            new_pair = [i,j-Md/2]
        else:
            continue
        
        bad_pairs.append(pair)
        good_pairs.append(new_pair)

    position_pairs = list(Combinations([t for t in range(dels)],Integer(2)))
    position_pairs.sort(key=lambda str: str[::-1])

    for tup in list_1:
        new_tup = deepcopy(tup)
        counts = [0 for i in range(dels)] # fix parity
        for posps in position_pairs:
            if [new_tup[posps[0]],new_tup[posps[1]]] in bad_pairs:
                good_pair = good_pairs[bad_pairs.index(
                    [new_tup[posps[0]],new_tup[posps[1]]])]
                new_tup[posps[0]] = good_pair[0]
                new_tup[posps[1]] = good_pair[1]

                counts[posps[0]] = counts[posps[0]] + 1
                counts[posps[1]] = counts[posps[1]] + 1
        
        if not is_weakly_decreasing_nonzero(counts): #dont know how to justify this/why it works
            for posps in reversed(position_pairs):
                if [new_tup[posps[0]],new_tup[posps[1]]] in bad_pairs:
                    good_pair = good_pairs[bad_pairs.index(
                        [new_tup[posps[0]],new_tup[posps[1]]])]
                    new_tup[posps[0]] = good_pair[0]
                    new_tup[posps[1]] = good_pair[1]
                    break
        
        while len(list(set(new_tup))) < dels:
            for posps in position_pairs:
                if [new_tup[posps[0]],new_tup[posps[1]]] in bad_pairs:
                    good_pair = good_pairs[bad_pairs.index(
                        [new_tup[posps[0]],new_tup[posps[1]]])]
                    new_tup[posps[0]] = good_pair[0]
                    new_tup[posps[1]] = good_pair[1]

                    
        list_2.append(new_tup)
        
    # print(bad_pairs)
    return(list_1,list_2)


m,k,dels = 4,0,3
a,b = CFI(K[m],0),CFI(K[m],1)
n = a.nrows()
switch = 0
Fd = n/m
Od = 2*(m-1)
Md = Fd-Od


# bij = complete_CFI_explicit_bijection_2(m,dels)
# for p in range(len(bij[0])):
#     if bij[0][p][0] == 10 and not list_del_tester(a,b,bij[0][p],bij[1][p],k):
#         print(bij[0][p],bij[1][p])

for i in range(n):
    print("round",i)
    for j in range(n):
        for kp in range(n):
            if (not list_del_tester(a,b,[Fd,kp,i],[Fd,kp,i],0)) and i not in [0,1,6,7,8,9]:
                if list_del_tester(a,b,[Fd,kp,i],[Fd,kp,j],0):
                    print(i,kp,j)

# print(list_del_tester(a,b,[10,20,j],[10,39,i],0))

# bijection = complete_CFI_explicit_bijection_2(m,dels)
# for i in range(len(bijection[0])):
#     if True:##practicalities
#         if not list_del_tester(a,b,bijection[0][i],bijection[1][i],k) or bijection[0][i][0]!=bijection[1][i][0]:
#             print(bijection[0][i],bijection[1][i],"is the problem")
#             switch = 1
#             break
#         else:
#             print(str(i)+"/"+str(len(bijection[0])-1),bijection[0][i],bijection[1][i])
#         if switch == 1:
#             break
#     else:####
#         print(str(i)+"/"+str(len(bijection[0])-1))