from copy import copy
from itertools import cycle


def is_k_periodic(lst, k):
    if len(lst) < k // 2:  # we want the returned part to repeat at least twice... otherwise every list is periodic (1 period of its full self)
        return False

    return all(x == y for x, y in zip(lst, cycle(lst[:k])))


def is_periodic(lst):
    for k in range(1, (len(lst) // 2) + 1):
        if is_k_periodic(lst, k):
            return tuple(lst[:k])
    return None

def is_sum_free(nums):
    for i in range(len(nums)):
        for j in range(i,len(nums)):
            k = nums[i]+nums[j]
            if k in nums:
                return(False)
    return(True)

def whole_cameron(nums, lim):
    if not is_sum_free(nums):
        raise TypeError("not sum-free")

    sumlist = copy(nums)
    old_len = 0
    counter = 1

    while counter < max(nums):
        for i in range(len(nums)):
            for j in range(i,len(nums)):
                k = nums[i]+nums[j]
                if k not in sumlist:
                    sumlist.append(k)
        
        if is_sum_free([counter]+nums) and counter not in nums:
            nums.append(counter)
        
        counter = counter + 1

    nums.sort()

    while counter < lim:
        if counter not in sumlist:
            for i in nums:
                if i + counter not in sumlist:
                    sumlist.append(i + counter)
            nums.append(counter)

        counter = counter + 1
    
    diffs = []
    for i in range(len(nums)-1):
        diffs.append(nums[i+1]-nums[i])

    print(nums)
    # print(diffs)
    print(is_periodic(diffs))
    return(diffs)

def follow_cameron(inums, lim):

    nums = copy(inums)
    if not is_sum_free(nums):
        raise TypeError("not sum-free")

    initial_len = len(nums)
    sumlist = copy(nums)
    counter = max(nums)

    while counter < lim:
        if counter not in sumlist:
            for i in nums:
                if i + counter not in sumlist:
                    sumlist.append(i + counter)
            nums.append(counter)

        counter = counter + 1
    
    for i in range(initial_len):
        del(nums[0])

    diffs = []
    for i in range(len(nums)-1):
        diffs.append(nums[i+1]-nums[i])

    print(inums,nums)
    # print(diffs)
    print(is_periodic(diffs))
    return(diffs)

# whole_cameron([1,7,9,12],100)
whole_cameron([1,7,9,12,32,51],200)
whole_cameron([1,12,32,51],200)

# follow_cameron([1,7,9,12],200)
follow_cameron([1,7,9,12,32,51],1000)
follow_cameron([1,12,32,51],1000)
is_periodic([2, 7, 2, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 2, 7, 2, 2, 7, 2, 2, 9, 2, 9, 2, 2])
