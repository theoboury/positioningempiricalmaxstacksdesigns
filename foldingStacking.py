# positioningempiricalmaxstacksdesigns
# Copyright (C) 2026 THEO BOURY 

#This file contains the function to launch folding algorithms. 
# #The energy models supported are unitary or stacking. Weighted can be added by modifying "Energy" function.
#Considered BPs can be G-U/U-G or not.
#Method is based on dynamic programming. Every time two backtracks are provided:
# - one leading deterministically to an optimal solution 
# - one leading to all solutions at Delta from the optimal.
#For now, theta (min between two extremities of BP) = 0 and is not reported in the function.
#A Varna interface is used to draw an optimal sequence and structure.

import math
        
def isValid(seq, i, j, BPconsidered="Nussinov"):
    """Input: - A Sequence seq of letters in {A, U, C, G} of size n
              - A BP (i, j) with i < n and j < n
              - A model of BP considered BPconsidered
       Output: If the BP (i, j) is feasible over seq or not from the "letters" constraints"""
    if BPconsidered == "Nussinov":
        S = [("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")]
    elif BPconsidered == "Watson":
        S = [("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")]
    else:
        raise ValueError("Not a valid model for BPs")
    return ((seq[i], seq[j]) in S)


def Energy(seq, i, j, model="Unitary", BPconsidered="Nussinov"):
    """Input: - A Sequence seq of letters in {A, U, C, G} of size n
              - A BP (i, j) with i < n and j < n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: The energy of BP (i,j) in the specified model, currently only unitary energetic contribution is supported"""
    if model == "Unitary":
        if isValid(seq, i, j, BPconsidered=BPconsidered):
            return -1
        return 0
    #Add models if wanted


def FillMatStacking2(seq, BPconsidered="Nussinov"):
    """Input: - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: The dynamic programming tables for folding over seq in a stacking model"""
    #No need to specify the model here, we need "unitary" contirbutions for each BP "stacked".
    n = len(seq)
    
    #W[i][j] is optimal energy for all structures possible between i and j
    E = [[0 for i in range(n)] for j in range(n)]
    
    #S[i][j] is optimal energy for all structures possible between i and j such that (i - 1,j + 1) is a valid BP. 
    S = [[0 for i in range(n)] for j in range(n)]

            
    for i in range(n - 1,-1, -1):
        for j in range(i + 1, n):
            #i base is unpaired
            case_unpaired = E[i + 1][j] 
            
            #(i, j) placed over (i+1, j-1) with possibility of stacking with.
            #if i != 0 and j !=n - 1:
            case_embraceS=math.inf
            case_embraceE=math.inf
            if isValid(seq, i, j, BPconsidered=BPconsidered):
                case_embraceS = S[i +1][j - 1] + Energy(seq, i, j, BPconsidered=BPconsidered) 
                case_embraceE = S[i +1][j - 1]
        
            #(i,k) forms a valid "pair" that split in a subinstance (i + 1, k - 1) and an exterior instance (k + 1, j).
            case_split = case_unpaired
            for k in range(i + 1, j):
                if isValid(seq, i, k, BPconsidered=BPconsidered):
                    case_split = min(case_split, S[i + 1][k - 1] + E[k + 1][j])
            E[i][j] = min(case_split, case_embraceE)
            S[i][j] = min(case_split, case_embraceS)
    return E, S


def DeltaBackTrackE2(sigma, Slist, delta, E, S, seq, model="Unitary", BPconsidered="Nussinov"):
    """Input: - Set of regions being considered sigma
              - Partial secondary structure asssigned to this point Slist
              - delta, max distance from the optimal allowed
              - The dynamic programming tables computed E and S
              - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: All optimal secondary structure over seq at distance delta from optimal"""
    
    if sigma == []:
        return [Slist]
    (i,j) = sigma.pop()

    if j<= i:
        return DeltaBackTrackE2(sigma, Slist, delta, E, S, seq, model=model, BPconsidered=BPconsidered)
    
    resu = []
    delt= S[i + 1][j - 1] - E[i][j]

    if delta - delt>= 0 and isValid(seq, i, j, BPconsidered=BPconsidered):

        newSlist = [(i, j)] + Slist
        newsigma = sigma + [(i + 1, j - 1)]
        resu += DeltaBackTrackS2(newsigma, newSlist, delta - delt, E, S, seq, model=model, BPconsidered=BPconsidered)
           
    delt= E[i + 1][j] - E[i][j]
    if delta - delt>= 0:
        newsigma = sigma + [(i + 1, j)]
        resu += DeltaBackTrackE2(newsigma, Slist, delta - delt, E, S, seq, model=model, BPconsidered=BPconsidered)
    
    for k in range(i + 1, j):
        if isValid(seq, i, k, BPconsidered=BPconsidered):
            delt =  S[i + 1][k -1] + E[k + 1][j] - E[i][j]
            if delta - delt>= 0:
                newsigma = sigma + [(k + 1, j), (i + 1, k - 1)]
                newSlist = [(i, k)] + Slist
                resu += DeltaBackTrackS2(newsigma, newSlist, delta - delt, E, S, seq, model=model, BPconsidered=BPconsidered)

    return resu

def DeltaBackTrackS2(sigma, Slist, delta, E, S, seq, model="Unitary", BPconsidered="Nussinov"):
    """Input: - Set of regions being considered sigma
              - Partial secondary structure asssigned to this point Slist
              - delta, max distance from the optimal allowed
              - The dynamic programming tables computed E and S
              - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: Aa optimal secondary structures over seq at delta from optimal, considering your are encapsulate by a  BP"""

    if sigma == []:
        return [Slist]
    (i,j) = sigma.pop()
    if j<= i:
        return DeltaBackTrackE2(sigma, Slist, delta, E, S, seq, model=model, BPconsidered=BPconsidered)
    
    resu = []
    delt= S[i + 1][j - 1] + Energy(seq, i, j, BPconsidered=BPconsidered) - S[i][j]
    if delta - delt>= 0 and isValid(seq, i, j, BPconsidered=BPconsidered):
        newSlist = [(i, j)] + Slist
        newsigma = sigma + [(i + 1, j - 1)]
        resu += DeltaBackTrackS2(newsigma, newSlist, delta - delt, E, S, seq, model=model, BPconsidered=BPconsidered)
           
    delt= E[i + 1][j] - S[i][j]
    if delta - delt>= 0:
        newsigma = sigma + [(i + 1, j)]
        resu += DeltaBackTrackE2(newsigma, Slist, delta - delt, E, S, seq, model=model, BPconsidered=BPconsidered)
    
    for k in range(i + 1, j):
        if isValid(seq, i, k, BPconsidered=BPconsidered):
            delt =  S[i + 1][k -1] + E[k + 1][j] - S[i][j]
            if delta - delt>= 0:
                newsigma = sigma + [(k + 1, j), (i + 1, k - 1)]
                newSlist = [(i, k)] + Slist
                resu += DeltaBackTrackS2(newsigma, newSlist, delta - delt, E, S, seq, model=model, BPconsidered=BPconsidered)
    
    return resu


def delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = "png", name_file = "", show=False):
    """Input: - A Sequence seq of letters in {A, U, C, G} of size n
             - A energy model to work with
             - A model of BP considered BPconsidered
             - delta, max distance from the optimal allowed
             - debug, print for debug
             - output_format is the type of file that we want to save, specify None for no save.
             - name_file is the name of the output Varna file if specified.
             - show is set to False if we want to observe the output directly or not (interactively).
      Output: Build a MFE structure, print it and draw it  in the Unitary model"""
    E, S = FillMatStacking2(seq, BPconsidered=BPconsidered)
    if debug:
        print(E)
        print(S)
    Slist = DeltaBackTrackE2([(0, len(seq) - 1)], [], delta, E, S, seq, model=model, BPconsidered=BPconsidered)
    for St in Slist:
        if debug:
            print(St)
    return Slist
