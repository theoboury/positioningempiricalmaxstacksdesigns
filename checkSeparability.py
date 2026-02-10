import random
from math import exp
def ssparse(seq):
    """
    Input:
        * seq, a secondary structure as a well-parenthesized string
    Output:
        * A list with -1 for an unpaired position and the index of its correspondent for a paired position
    """
    res = [-1 for c in seq]
    p = []
    for i,c in enumerate(seq):
        if c=='(':
            p.append(i)
        elif c== ')':
            j = p.pop()
            (res[i],res[j]) = (j,i)
    return res

def dbn_to_tree(dbn,i=None,j=None):
    """
    Input:
       * dbn, a list with -1 for an unpaired position and the index of its correspondent for a paired position
       * i, j, ints that corresponds to the encapsulating base pair that starts the tree (if exists) 
    Output:
       * A boolean if the node v a leaf or not 
    """
    if i is None:
        (i,j)=(0,len(dbn)-1)
        return ("root",dbn_to_tree(dbn,i,j))
    res = []
    while i<=j:
        if dbn[i] == -1:
            res.append(((i,i),[]))
            i+=1
        else:
            (ii,jj) = (i,dbn[i])
            res.append(((ii,jj),dbn_to_tree(dbn,ii+1,jj-1)))
            i = dbn[i]+1
    return res   

def is_leaf(v):
    """
    Input:
       * v, a secondary structure tree
    Output:
       * A boolean if the node v[0] is a leaf or not 
    """
    if v[0] == "root":
        return False
    (a, b) = v[0]
    if a == b:
        return True 
    return False

def children(v):
    """
    Input:
       * v, a secondary structure tree
    Output:
       * List of children of node v[0]
    """
    return v[1]

def filter(v):
    """
    Input:
       * v, a secondary structure tree
    Output:
       * A boolean if the tree contains no m3o and m5 or not 
    """
    BP_children = [vv[0] for vv in v[1] if not is_leaf(vv)]
    BP_leaves = [vv[0] for vv in v[1] if is_leaf(vv)]
    val = 3
    val2 = 1
    if v[0] == "root":
        val = 4
        val2 = 2
    if len(BP_children) > val or (len(BP_leaves) > 0 and len(BP_children) > val2):
        return False
    resu = True
    for vv in v[1]:
        resu = resu and filter(vv)
    return resu

children_colors_from_parent = {
    'AU' : ['AU','GC','CG'],
    'UA' : ['UA','GC','CG'],
    'GC' : ['UA','GC','AU'],
    'CG' : ['UA','CG','AU'],
    'R' : ['GC', 'UA','CG','AU']
}

def subset(li1, li2):
    """
    Input:
        * li1 and li2, two lists
    Output:
        * If elements of li1 (not necessarily distinct in li1) are a subsets of li2
    """
    sli1 = set(li1)
    bisli1 = list(sli1)
    bisli1.sort()
    if li1 == bisli1 and sli1.issubset(set(li2)):
        return True
    return False

def isProper(v, seq):
    """
    Input:
        * v, a tree 
        * seq, a sequence 
    Output:
        * A boolean that says if the coloring associated with seq is proper
    """
    if is_leaf(v):
        return (seq[v[0][0]] == "A")
    else:
        if v[0] == "root":
            tolerated_colors = children_colors_from_parent['R']
        else:
            c = seq[v[0][0]] + seq[v[0][1]]
            if c not in ["AU", "UA", "GC", "CG"]:
                return False
            tolerated_colors = children_colors_from_parent[c]
        BP_children = [vv[0] for vv in v[1] if not is_leaf(vv)]
        colors_BP = [seq[name[0]] + seq[name[1]] for name in BP_children]
        colors_BP.sort()
        if not subset(colors_BP, tolerated_colors):
            return False
        else:
            resu = True
            for w in children(v):
                resu = resu and isProper(w, seq)
            return resu


def isSeparable(t, seq, minmodulo = -1):
    """
    Input:
        * t, a tree
        * seq, a sequence that corresponds to a proper coloring
    Output:
        * A boolean that says if the coloring associated with seq is separable
          This coloring is supposed to be proper
    """
    def leaf_to_A(v):
        if is_leaf(v):
            if seq[v[0][0]] != "A":
                return False
            else:
                return True
        else:
            resu = True
            for w in children(v):
               resu = resu and leaf_to_A(w) 
            return resu
    def grey_and_leaf_level(v, current):
        if is_leaf(v):
            return [("l", current)]
        else:
            resu = []
            next_level = current
            if v[0] != "root":
                c = seq[v[0][0]] + seq[v[0][1]]
                if c == "AU" or c=="UA":
                    resu = [("g", current)]
                elif c == "GC":
                    next_level = current + 1  
                elif c == "CG":
                    next_level = current - 1  
            for w in children(v):
               resu = resu + grey_and_leaf_level(w, next_level)    
            return resu
    if not leaf_to_A(t):
        if minmodulo == -1:
            return False
        return False, -1
    li = grey_and_leaf_level(t, 0)
    LV_grey = [j for (i,j) in li if i == "g"]
    LV_leaf = [j for (i,j) in li if i == "l"]
    inter = [j for j in LV_leaf if j in LV_grey]
    m = -1
    for i in range(2, minmodulo + 1):
        LV_grey_mod = [j%i for j in LV_grey]
        LV_leaf_mod = [j%i for j in LV_leaf]
        inter_mod = [j for j in LV_leaf_mod if j in LV_grey_mod]
        if (inter_mod == []):
            m = i
            break
    if minmodulo == -1:
        return (inter == [])
    else:
        return (inter == []), m


def fullSeparable(seq, ss):
    t = dbn_to_tree(ssparse(ss))
    if filter(t):
        if isProper(t, seq):
            if isSeparable(t, seq):
                return True
    return False

def delta(c):
    """
    Input:
        * c, assignment of the base pairs (or information that it is the root)
    Output:
        * The difference of levels after the node assigned c
    """
    if c=="GC":
        return 1
    elif c=="CG":
        return -1
    elif c=="AU":
        return 0
    elif c=="UA":
        return 0
    elif c == "R":
        return 0
    raise("Kernel panic!")
    



def get_assignments(l,c,i=0):
    """
    Input:
        * l, a (partial) list of an ordered list of assignments for the leaves and BP children
        * c, assignment of the parent base pairs (or information that the parent is the root)
        * i, the index of the node considered to assign
    Output:
        * The complete list of an ordered list of assignments for the leaves and BP children
    """
    if i >= len(l):
        return [[]]
    else:
        v = l[i]
        if is_leaf(v):
            if c not in ["AU","UA"]:
                return [['A'] + lp for lp in get_assignments(l,c,i+1) if ('AU' not in lp) and ('UA' not in lp)]
            else:
                return []
        
        else:
            res = []
            for cv in children_colors_from_parent[c]:
                for a in get_assignments(l,c,i+1):
                    if cv not in a:
                        if (cv not in ["AU","UA"]) or ('A' not in a):
                            res.append([cv] + a)
            return res


    
    
def num_design(v, c, current_level_mod, m, leaves_levels_mod, cache, GCweight=None):
    """
    Input:
        * v, a secondary structure tree
        * c, the assignement for node v[0]
        * current_level_mod, the current level modulo m
        * m, the modulo considered for m-separability
        * leaves_levels_mod, the list of levels specific to the nodes
        * cache, the partial solutions for subtrees
        * GCweight, optional, the weight to attribute to GC base pairs
    Output:
        * The number of designs for the current tree under these parameters
    """
    state = (v[0], c, current_level_mod,m,leaves_levels_mod)
    if state not in cache:
        if is_leaf(v):
            if current_level_mod in leaves_levels_mod:
                cache[state] = 1
            else:
                cache[state] = 0
        elif (c=="AU" or c=="UA") and (current_level_mod in leaves_levels_mod):
            cache[state] = 0
        else:
            acc = 0
            next_level_mod = (current_level_mod + delta(c)) % m
            for assignment in get_assignments(children(v),c):
                prod = 1
                if GCweight is not None and (c in ["GC", "CG"]):
                    prod = exp(GCweight)
                for i,w in enumerate(children(v)):
                    prod *= num_design(w, assignment[i], next_level_mod, m,leaves_levels_mod, cache, GCweight=GCweight)
                acc += prod
            cache[state] = acc
    return cache[state]


def stochastic_backtrack(v, c, current_level_mod, m, leaves_levels_mod, cache, GCweight=None):
    """
    Input:
        * v, a secondary structure tree
        * c, the assignement for node v[0]
        * current_level_mod, the current level modulo m
        * m, the modulo considered for m-separability
        * leaves_levels_mod, the list of levels specific to the nodes
        * cache, the partial solutions for subtrees
        * GCweight, optional, the weight to attribute to GC base pairs
    Output:
        * A design uniformly sampled for the current tree under these parameters
    """
    state = (v[0], c, current_level_mod, m, leaves_levels_mod)
    if is_leaf(v):
        return c
    else:
        acc = 0
        next_level_mod = (current_level_mod + delta(c)) % m
        r = random.random()*cache[state]
        for assignment in get_assignments(children(v),c):
            prod = 1
            if GCweight is not None and (c in ["GC", "CG"]):
                prod = exp(GCweight)
            for i,w in enumerate(children(v)):
                prod *= num_design(w, assignment[i], next_level_mod, m,leaves_levels_mod, cache, GCweight=GCweight)
            r -= prod
            if r<0:
                res = ""
                for i,w in enumerate(children(v)):
                    loc = stochastic_backtrack(w, assignment[i], next_level_mod, m,leaves_levels_mod, cache, GCweight=GCweight)
                    if loc is None:
                        return None
                    res += loc
                if c == "R":
                    return res
                else:
                    return c[0]+res+c[1]


def part(i):
    """
    Input:
        * i, an integer
    Output:
        * A list of all distinct sets of i + 1 elements (as lists)
    """
    if i == -1:
        return [[]]
    else:
        li = part(i - 1)
        return li + [elem+[i] for elem in li]


    






def first_modulo_separable(t, modulolimit=4):
    for m in range(2, modulolimit + 1):
        for leaf_level in part(m - 1):
            cache = {}
            if num_design(t, "R", 0, m, tuple(leaf_level), cache, GCweight=None) > 0:
                seq = stochastic_backtrack(t, "R", 0, m, tuple(leaf_level), cache)
                return seq
            