from copy import copy

def circulant(nrow, ncol, vec):
    if len(vec) != ncol:
        raise TypeError("Size mismatch")
    M = []
    for i in range(nrow):
        M.append(copy(vec))
        tmp = vec[-1]
        vec.pop()
        vec = [tmp] + vec
    return(M)

print(circulant(3,4,[1,2,3,4]))

list_3_names = [26,28,"35 18",36,40,45,50]
for i in range(len(list_3_names)):
    list_3_names[i]

p = [3]
7<min(p)

p = [[3], [4], [3], [3]]
p.sort()
print(p)

foot = [1,1,2,3]
colour = []
uniques = list(set(foot))
for val in uniques:
    colour.append(foot.count(val))
colour.sort()
print(colour)

for i in range(10, -1, -1):
    print(i)