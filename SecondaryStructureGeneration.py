# positioningempiricalmaxstacksdesigns
# Copyright (C) 2026 THEO BOURY 

import random 
UNPAIRED_WEIGHT = 1

#Below are some functions for the random generation of secondary structures
def sscount(n,count,count_stacked,theta,min_helix):
    """
    Input:
        * n, size of the secondary structures
        * count, a dictionary of the number of secondary structures of a given size. The size serves as an index
        * count_stacked, same as count but count only secondary structures finishing by a base pair
        * theta, the minimal number of unpaired positions between the extremities of a base pair
        * min_helix, minimal size allowed for the helices
    Output:
        * An int, number of secondary structures of size n
    """
    if n not in count:
        if n<=0:
            count[n] = 1.
        else:
            count[n] = UNPAIRED_WEIGHT*sscount(n-1,count,count_stacked,theta,min_helix)
            if n >= theta+2*min_helix:
                count[n] += sscount_stacked(n-2*min_helix,count,count_stacked,theta,min_helix)
            for i in range(theta+2*min_helix,n):
                count[n] += sscount_stacked(i-2*min_helix,count,count_stacked,theta,min_helix)*sscount(n-i,count,count_stacked,theta,min_helix)
    return count[n]


def sscount_stacked(n,count,count_stacked,theta,min_helix):
    """
    Input:
        * n, size of the secondary structures
        * count, a dictionary of the number of secondary structures of a given size. The size serves as an index
        * count_stacked, same as count but count only secondary structures finishing by a base pair
        * theta, the minimal number of unpaired positions between the extremities of a base pair
        * min_helix, minimal size allowed for the helices
    Output:
        * An int, number of secondary structures that finish by a base pair of size n
    """
    if n not in count_stacked:
        if n<=0:
            count_stacked[n] = 1.
        else:
            count_stacked[n] = UNPAIRED_WEIGHT*sscount(n-1,count,count_stacked,theta,min_helix)
            if n >= theta+2:
                count_stacked[n] += sscount_stacked(n-2,count,count_stacked,theta,min_helix)
            for i in range(theta+2*min_helix,n):
                count_stacked[n] += sscount_stacked(i-2*min_helix,count,count_stacked,theta,min_helix)*sscount(n-i,count,count_stacked,theta,min_helix)
    return count_stacked[n]


def ssrandom(n,count,count_stacked,theta,min_helix):
    """
    Input:
        * n, size of the secondary structures
        * count, a dictionary of the number of secondary structures of a given size. The size serves as an index
        * count_stacked, same as count but count only secondary structures finishing by a base pair
        * theta, the minimal number of unpaired positions between the extremities of a base pair
        * min_helix, minimal size allowed for the helices
    Output:
        * A uniformly sampled secondary structure
    """
    if n<=0:
        return ''
    else:
        r = random.random()*sscount(n,count,count_stacked,theta,min_helix)
        r -= UNPAIRED_WEIGHT*sscount(n-1,count,count_stacked,theta,min_helix)
        if r<0:
            return '.'+ssrandom(n-1,count,count_stacked,theta,min_helix)
        if n >= theta+2*min_helix:
            r -= sscount_stacked(n-2*min_helix,count,count_stacked,theta,min_helix)
            if r<0:
                return '('*min_helix+ssrandom_stacked(n-2*min_helix,count,count_stacked,theta,min_helix)+')'*min_helix
        for i in range(theta+2*min_helix,n):
            r -= sscount_stacked(i-2*min_helix,count,count_stacked,theta,min_helix)*sscount(n-i,count,count_stacked,theta,min_helix)
            if r<0:
                return '('*min_helix+ssrandom_stacked(i-2*min_helix,count,count_stacked,theta,min_helix)+')'*min_helix+ssrandom(n-i,count,count_stacked,theta,min_helix)


def ssrandom_stacked(n,count,count_stacked,theta,min_helix):
    """   
    Input:
        * n, size of the secondary structures
        * count, a dictionary of the number of secondary structures of a given size. The size serves as an index
        * count_stacked, same as count but count only secondary structures finishing by a base pair
        * theta, the minimal number of unpaired positions between the extremities of a base pair
        * min_helix, minimal size allowed for the helices
    Output:
        * A uniformly sampled secondary structure that finish by a base pair
    """ 
    if n<=0:
        return ''
    else:
        r = random.random()*sscount_stacked(n,count,count_stacked,theta,min_helix)
        r -= UNPAIRED_WEIGHT*sscount(n-1,count,count_stacked,theta,min_helix)
        if r<0:
            return '.'+ssrandom(n-1,count,count_stacked,theta,min_helix)
        if n >= theta+2:
            r -= sscount_stacked(n-2,count,count_stacked,theta,min_helix)
            if r<0:
                return '('+ssrandom_stacked(n-2,count,count_stacked,theta,min_helix)+')'
        for i in range(theta+2*min_helix,n):
            r -= sscount_stacked(i-2*min_helix,count,count_stacked,theta,min_helix)*sscount(n-i,count,count_stacked,theta,min_helix)
            if r<0:
                return '('*min_helix+ssrandom_stacked(i-2*min_helix,count,count_stacked,theta,min_helix)+')'*min_helix+ssrandom(n-i,count,count_stacked,theta,min_helix)
