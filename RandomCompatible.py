import random
def children(v):
    """
    Input:
       * v, a secondary structure tree
    Output:
       * List of children of node v[0]
    """
    return v[1]

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

def choose_random_seq(t, withA=True):
    """
    Input:
        * t, a tree 
        * withA, a boolean to say if we place only A on the unpaired positions or not
    Output:
        * A random sequence compatible with the tree t
    """
    def seq_from_tree(v):
        if is_leaf(v):
            if withA:
                return "A"
            else:
                return random.choice(["A", "G", "C", "U"])
        else:
            res = ""
            BP = random.choice(["AU","UA", "GC", "CG"])
            for w in children(v):
                res= res + seq_from_tree(w)
            if v[0] == "root":
                return res
            else:
                return BP[0] + res + BP[1]
    return seq_from_tree(t)
