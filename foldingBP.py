# positioningempiricalmaxstacksdesigns
# Copyright (C) 2026 THEO BOURY 

#This file contains the function to launch folding algorithms. 
# #The energy model supported is unitary. Weighted can be added by modifying the "Energy" function.
#Considered BPs can contain G-U/U-G or not.
#For now, theta (min between two extremities of BP) = 0 and is not reported in the function.


def isValid(seq, i, j, BPconsidered="Nussinov"):
    """
    Input: 
        * A Sequence seq of letters in {A, U, C, G} of size n
        * A BP (i, j) with i < n and j < n
        * A model of BP considered BPconsidered
    Output: 
        * If the BP (i, j) is feasible over seq or not from the "letters" constraints
    """
    if BPconsidered == "Nussinov":
        S = [("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")]
    elif BPconsidered == "Watson":
        S = [("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")]
    else:
        raise ValueError("Not a valid model for BPs")
    return ((seq[i], seq[j]) in S)


def Energy(seq, i, j, model="Unitary", BPconsidered="Nussinov"):
    """
    Input: 
        * A Sequence seq of letters in {A, U, C, G} of size n
        * A BP (i, j) with i < n and j < n
        * An energy model to work with
        * A model of BP considered BPconsidered
    Output: 
        * The energy of BP (i,j) in the specified model, currently only unitary energetic contribution is supported
    """
    if model == "Unitary":
        if isValid(seq, i, j, BPconsidered=BPconsidered):
            return -1
        return 0
    #Add models if wanted


def FillMatUnitary(seq, model="Unitary", BPconsidered="Nussinov"):
    """
    Input: 
        * A sequence seq of nucleotides {A, U, C, G} of size n
        * An energy model to work with
        * A model of BP considered BPconsidered
    Output: 
        * The dynamic programming table for folding over seq in a unitary or weighted model
    """
    n = len(seq)
    M = [[0 for i in range(n)] for j in range(n)]
    N = [[1 for i in range(n)] for j in range(n)]       
    for i in range(n - 1,-1, -1):
        for j in range(i + 1, n):
            
            #i + 1 base is unpaired
            case_unpaired = M[i + 1][j] 
            
            #(i,j) is a "pair". No need to distinguish cases depending if the pair is valid or not, not the case with stackings
            case_embrace = M[i + 1][j - 1] + Energy(seq, i, j, model="Unitary", BPconsidered="Nussinov") 
            
            #(i,k) forms a "pair" (valid or not) that splits in a subinstance (i + 1, k - 1) and an exterior instance (k + 1, j).
            case_split = min(case_unpaired, case_embrace) 
            for k in range(i + 1, j):
                case_split = min(case_split, M[i + 1][k - 1] + M[k + 1][j] + Energy(seq, i, k))
            M[i][j] = case_split

            N[i][j] = 0
            if M[i][j] == case_unpaired:
                N[i][j] += N[i + 1][j] 
            if M[i][j] == case_embrace:
                N[i][j] += N[i + 1][j - 1] 
            for k in range(i + 1, j):
                if M[i][j] == M[i + 1][k - 1] + M[k + 1][j] + Energy(seq, i, k):
                   N[i][j] += N[i + 1][k - 1] + N[k + 1][j]
    return M, N[0][n - 1]
     

def DeltaBackTrackUnitary(sigma, Slist, delta, M, seq, model="Unitary", BPconsidered="Nussinov"):
    """
    Input: 
        * Set of regions being considered sigma
        * Partial secondary structure assigned to this point Slist
        * delta, max distance from the optimal allowed
        * The dynamic programming table computed M
        * A sequence seq of nucleotides {A, U, C, G} of size n
        * An energy model to work with
        * A model of BP considered BPconsidered
    Output: 
        * All optimal secondary structures over seq at delta from the optimal
    """
       
    if sigma == []:
        return [Slist]
    
    (i,j) = sigma.pop()
    if j<= i:
        return DeltaBackTrackUnitary(sigma, Slist, delta, M, seq, model=model, BPconsidered=BPconsidered)
    
    resu = []
    
    delt= M[i + 1][j] - M[i][j]
    if delta - delt>= 0:
        newsigma = [(i + 1, j)] + sigma
        resu += DeltaBackTrackUnitary(newsigma, Slist, delta - delt, M, seq, model=model, BPconsidered=BPconsidered)
    
    delt= M[i + 1][j - 1] + Energy(seq, i, j, model=model, BPconsidered=BPconsidered) - M[i][j]
    if delta - delt>= 0 and isValid(seq, i, j, BPconsidered=BPconsidered):
        newsigma = [(i + 1, j - 1)] + sigma
        newSlist = [(i, j)] + Slist
        resu += DeltaBackTrackUnitary(newsigma, newSlist, delta - delt, M, seq, model=model, BPconsidered=BPconsidered)
    

    for k in range(i + 1, j):
        delt= M[i + 1][k -1] + M[k + 1][j] + Energy(seq, i, k, model=model, BPconsidered=BPconsidered) - M[i][j]
        if delta - delt>= 0 and isValid(seq, i, k, BPconsidered=BPconsidered):
            newsigma = [(i + 1, k - 1), (k + 1, j)] + sigma
            newSlist = [(i, k)] + Slist
            resu += DeltaBackTrackUnitary(newsigma, newSlist, delta - delt, M, seq, model=model, BPconsidered=BPconsidered)
    return resu


def OptimalBackTrackUnitary(i, j, M, seq, model="Unitary", BPconsidered="Nussinov"):
    """Input: - Extremities of  studied region[i, j]
              - The dynamic programming table computed M
              - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: An optimal secondary structure over seq"""
    if j<= i:
        return []
    if M[i][j] == M[i + 1][j]:
        return OptimalBackTrackUnitary(i + 1, j, M, seq, model=model, BPconsidered=BPconsidered)
    if M[i][j] == M[i + 1][j - 1] + Energy(seq, i, j, model=model, BPconsidered=BPconsidered): #WARNING: Here no need to distinguish case when the BP is valid or not, it will change with stackings
        addon = []
        if isValid(seq, i, j, BPconsidered=BPconsidered):
            addon = [(i, j)]
        return addon + OptimalBackTrackUnitary(i + 1, j - 1, M, seq, model=model, BPconsidered=BPconsidered)
    for k in range(i + 1, j):
        if M[i][j] == M[i + 1][k -1] + M[k + 1][j] + Energy(seq, i, k, model=model, BPconsidered=BPconsidered):
            addon = []
            if isValid(seq, i, k, BPconsidered=BPconsidered):
                addon = [(i, k)]
            addon += OptimalBackTrackUnitary(i + 1, k - 1, M, seq, model=model, BPconsidered=BPconsidered)
            addon += OptimalBackTrackUnitary(k + 1, j, M, seq, model=model, BPconsidered=BPconsidered)
            return addon
        

def main_unitary_only_one(seq, model="Unitary", BPconsidered="Nussinov"):
    """
    Input: 
        * A Sequence seq of letters in {A, U, C, G} of size n
        * A energy model to work with
        * A model of BP considered BPconsidered
        * delta, max distance from the optimal allowed
    Output: 
        * One structure that maximize the number of base pairs and if there are multiple structures
    """
    M, nb = FillMatUnitary(seq, model=model, BPconsidered=BPconsidered)
    S = OptimalBackTrackUnitary(0, len(seq) - 1, M, seq, model=model, BPconsidered=BPconsidered)
    return S, nb

