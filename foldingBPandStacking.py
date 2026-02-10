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

import varnaapi
import math
from random import seed, choice

def convert(seq, listS, BPconsidered="Nussinov", output_format = "png", name_file = "", show=False, seeder=42, full = 0):
    """Input: - A sequence seq of nucleotides {A, U, C, G} of size n
              - A list of list of BPs S over seq
              - A model of BP considered BPconsidered
              - output_format is the type of file that we want to save, specify None for no save.
              - name_file is the name of the output Varna file if specified.
              - show is set to False if we want to observe the output directly or not (interactively).
              - a seed for random generation
       Output: If it exists, returns the well parenthesized structure for S over seq"""
    #seed(seeder)
    n = len(seq)
    #print(seq)
    list_resu = []
    for S in listS:
        #print(S)
        resu = ["." for i in range(n)]
        for (i,j) in S:
            if not isValid(seq, i , j, BPconsidered=BPconsidered):
                return None
            resu[i] = "("
            resu[j] = ")"
        resu= "".join(resu)
        if resu != "."*n:
            list_resu.append(resu)
        #print(resu)
        #print("\n")
    #print(listS)
    if S != []:
        S = choice(listS)
        v  = varnaapi.Structure(seq, S)
        if output_format:
            if name_file == "":
                name = "SoverSeq" + "." + output_format
            else:
                name = name_file + "." + output_format
            v.savefig(name)
        if show:
            v.show()
    #print("hola")
    if full:
        return list_resu
    return len(listS)
        
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


def FillMatUnitary(seq, model="Unitary", BPconsidered="Nussinov"):
    """Input: - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: The dynamic programming table for folding over seq in a unitary or weigthed model"""
    n = len(seq)
    M = [[0 for i in range(n)] for j in range(n)]
    #for i in range(n): #No use here for this model, but can be useful depending on the initialization of the assignations
    #    for j in range(i, n): #min(i, n - 1) + 1 when theta <> 0
    #        M[i][j] = 0 
            
    for i in range(n - 1,-1, -1):
        for j in range(i + 1, n):
            
            #i + 1 base is unpaired
            case_unpaired = M[i + 1][j] 
            
            #(i,j) is a "pair". No need to distinguish cases depending if the pair is valid or not, not the case with stackings
            case_embrace = M[i + 1][j - 1] + Energy(seq, i, j, model="Unitary", BPconsidered="Nussinov") 
            
            #(i,k) forms a "pair" (valid or not) that split in a subinstance (i + 1, k - 1) and an exterior instance (k + 1, j).
            case_split = min(case_unpaired, case_embrace) 
            for k in range(i + 1, j):
                case_split = min(case_split, M[i + 1][k - 1] + M[k + 1][j] + Energy(seq, i, k))
            M[i][j] = case_split
            
    return M

def FillMatStacking(seq, BPconsidered="Nussinov"):
    """Input: - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: The dynamic programming tables for folding over seq in a stacking model"""
    #No need to specify the model here, we need "unitary" contirbutions for each BP "stacked".
    n = len(seq)
    
    #W[i][j] is optimal energy for all structures possible between i and j
    W = [[0 for i in range(n)] for j in range(n)]
    
    #S[i][j] is optimal energy for all structures possible between i and j such as (i,j) is a valid BP. 
    S = [[math.inf for i in range(n)] for j in range(n)]

            
    for i in range(n - 1,-1, -1):
        for j in range(i + 1, n):
            #S case
            if isValid(seq, i, j, BPconsidered=BPconsidered): #make sure that the S case is possible, kept infinite otherwise
                #(i, j) stacked over (i+1, j-1), case_unstack will always be prefered if the nergy is 0n, no need to case over validity of (i+1, j-1)
                case_stack = S[i +1][j - 1] + Energy(seq, i, j, BPconsidered=BPconsidered) 
        
                #(i, j) is a BP but does not had a "stacking" contribution
                case_unstack = W[i + 1][j - 1]
                S[i][j] = min(case_stack, case_unstack)
            
            #W case
            #i + 1 base is unpaired
            case_unpaired = W[i + 1][j]
            
            #(i, j) is a valid BP, use S for that
            caseS = S[i][j]
            #(i,k) forms a "pair" (valid or not) that split in a subinstance (i + 1, k - 1) and an exterior instance (k + 1, j).
            case_split = min(case_unpaired, caseS)
            for k in range(i + 1, j):
                if isValid(seq, i, k, BPconsidered=BPconsidered):
                    w = W[k + 1][j]
                    case_split = min(case_split, S[i + 1][k - 1] + w + Energy(seq, i, k)) #(i,k) contributes
                    case_split = min(case_split, W[i + 1][k - 1] + w) #(i,k) is not stacked, this case is always prefered if (i + 1, k - 1) non valid
            W[i][j] = case_split
    return W, S

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
        
        
def OptimalBackTrackW(i, j, W, S, seq, model="Unitary", BPconsidered="Nussinov"):
    #WARNING: You should prefer OptimalBackTrackE2 instead that is non ambiguous
    """Input: - Extremities of  studied region[i, j]
              - The dynamic programming tables computed W and S
              - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: An optimal secondary structure over seq"""
    if j<= i:
        return []
    if W[i][j] == S[i][j]:
        return OptimalBackTrackS(i, j, W, S, seq, model=model, BPconsidered=BPconsidered)
    if W[i][j] == W[i + 1][j]:
        return OptimalBackTrackW(i + 1, j, W, S, seq, model=model, BPconsidered=BPconsidered)

    for k in range(i + 1, j):
        if isValid(seq, i, k, BPconsidered=BPconsidered):
            addon = [(i, k)]
            if W[i][j] == S[i + 1][k -1] + W[k + 1][j] + Energy(seq, i, k, model=model, BPconsidered=BPconsidered): 
                addon += OptimalBackTrackS(i + 1, k - 1, W, S, seq, model=model, BPconsidered=BPconsidered)
                addon += OptimalBackTrackW(k + 1, j, W, S, seq, model=model, BPconsidered=BPconsidered)
                return addon
            elif W[i][j] == W[i + 1][k -1] + W[k + 1][j]:
                addon += OptimalBackTrackW(i + 1, k - 1, W, S, seq, model=model, BPconsidered=BPconsidered)
                addon += OptimalBackTrackW(k + 1, j, W, S, seq, model=model, BPconsidered=BPconsidered)
                return addon


def OptimalBackTrackS(i, j, W, S, seq, model="Unitary", BPconsidered="Nussinov"):
    #WARNING: You should prefer OptimalBackTrackS2 instead that is non ambiguous
    """Input: - Extremities of  studied region[i, j]
              - The dynamic programming tables computed W and S
              - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: An optimal secondary structure over seq with a BP between i and j"""
    if j<= i:
        return []
    addon = [(i, j)]
    print("addon", addon[0])
    if S[i][j] == S[i + 1][j - 1] + Energy(seq, i, j, model=model, BPconsidered=BPconsidered):
        return addon + OptimalBackTrackS(i + 1, j - 1, W, S, seq, model=model, BPconsidered=BPconsidered)
    if S[i][j] == W[i + 1][j - 1]:
        return addon + OptimalBackTrackW(i + 1, j - 1, W, S, seq, model=model, BPconsidered=BPconsidered)

def OptimalBackTrackE2(i, j, E, S, seq, model="Unitary", BPconsidered="Nussinov"):
    """Input: - Extremities of  studied region[i, j]
              - The dynamic programming tables computed E and S
              - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: An optimal secondary structure over seq"""
       
    if j<= i:
        return []
    
    if E[i][j] == E[i + 1][j]:
        return OptimalBackTrackE2(i + 1, j, E, S, seq, model=model, BPconsidered=BPconsidered)
    
    if E[i][j] == S[i + 1][j - 1] and isValid(seq, i, j, BPconsidered=BPconsidered):
        addon = [(i, j)]
        return addon + OptimalBackTrackS2(i + 1, j - 1, E, S, seq, model=model, BPconsidered=BPconsidered)
    
    for k in range(i + 1, j):
        if isValid(seq, i, k, BPconsidered=BPconsidered):
            addon = [(i, k)]
            if E[i][j] == S[i + 1][k -1] + E[k + 1][j]:
                addon += OptimalBackTrackS2(i + 1, k - 1, E, S, seq, model=model, BPconsidered=BPconsidered)
                addon += OptimalBackTrackE2(k + 1, j, E, S, seq, model=model, BPconsidered=BPconsidered)
                return addon

def OptimalBackTrackS2(i, j, E, S, seq, model="Unitary", BPconsidered="Nussinov"):
    """Input: - Extremities of  studied region[i, j]
              - The dynamic programming tables computed E and S
              - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: An optimal secondary structure over seq knowing that current structure is encapsulated between i - 1 and j + 1"""
    if j<= i:
        return []
    if S[i][j] == E[i + 1][j]:
        return OptimalBackTrackE2(i + 1, j, E, S, seq, model=model, BPconsidered=BPconsidered)
    if S[i][j] == S[i + 1][j - 1] + Energy(seq, i, j, BPconsidered=BPconsidered) and isValid(seq, i, j, BPconsidered=BPconsidered):
        addon = [(i, j)]
        return addon + OptimalBackTrackS2(i + 1, j - 1, E, S, seq, model=model, BPconsidered=BPconsidered)
    for k in range(i + 1, j):
        if isValid(seq, i, k, BPconsidered=BPconsidered):
            addon = [(i, k)]
            if S[i][j] == S[i + 1][k -1] + E[k + 1][j]:
                addon += OptimalBackTrackS2(i + 1, k - 1, E, S, seq, model=model, BPconsidered=BPconsidered)
                addon += OptimalBackTrackE2(k + 1, j, E, S, seq, model=model, BPconsidered=BPconsidered)
                return addon

def DeltaBackTrackUnitary(sigma, Slist, delta, M, seq, model="Unitary", BPconsidered="Nussinov"):
    """Input: - Set of regions being considered sigma
              - Partial secondary structure asssigned to this point Slist
              - delta, max distance from the optimal allowed
              - The dynamic programming table computed M
              - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: All optimal secondary structures over seq at delta from the optimal"""
       
    if sigma == []:
        return [Slist]
    
    (i,j) = sigma.pop()
    #print("i,j", i, j, "\n")
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

def DeltaBackTrackW(sigma, Slist, delta, W, S, seq, model="Unitary", BPconsidered="Nussinov"):
    #Should prefer DeltaBackTrackE2, that is non ambiguous
    """Input: - Set of regions being considered sigma
              - Partial secondary structure asssigned to this point Slist
              - delta, max distance from the optimal allowed
              - The dynamic programming tables computed W and S
              - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: All optimal secondary structures over seq at delta from the optimal"""

    if sigma == []:
        return [Slist]
    (i,j) = sigma.pop()
    if j<= i:
        return DeltaBackTrackW(sigma, Slist, delta, W, S, seq, model=model, BPconsidered=BPconsidered)
    
    resu = []
    delt= S[i][j] - W[i][j]
    if delta - delt>= 0:

        newsigma = sigma + [(i, j)]
        resu += DeltaBackTrackS(newsigma, Slist, delta - delt, W, S, seq, model=model, BPconsidered=BPconsidered)
           
    delt= W[i + 1][j] - W[i][j]
    if delta - delt>= 0:
        newsigma = sigma + [(i + 1, j)]
        resu += DeltaBackTrackW(newsigma, Slist, delta - delt, W, S, seq, model=model, BPconsidered=BPconsidered)
    
    for k in range(i + 1, j):
        if isValid(seq, i, k, BPconsidered=BPconsidered):
            delt1 =  S[i + 1][k -1] + W[k + 1][j] + Energy(seq, i, k, model=model, BPconsidered=BPconsidered) - W[i][j]
            delt2 = W[i + 1][k -1] + W[k + 1][j] - W[i][j]

            if delta - delt1>= 0 and isValid(seq, i + 1, k - 1, BPconsidered=BPconsidered):
                newsigma = sigma + [(k + 1, j), (i + 1, k - 1)]
                newSlist = [(i, k)] + Slist

                resu += DeltaBackTrackS(newsigma, newSlist, delta - delt1, W, S, seq, model=model, BPconsidered=BPconsidered)
            if delta - delt2>= 0:

                newsigma = sigma + [(i + 1, k - 1), (k + 1, j)]
                newSlist = [(i, k)] + Slist
                resu += DeltaBackTrackW(newsigma, newSlist, delta - delt2, W, S, seq, model=model, BPconsidered=BPconsidered)
    return resu

def DeltaBackTrackS(sigma, Slist, delta, W, S, seq, model="Unitary", BPconsidered="Nussinov"):
    """Input: - Set of regions being considered sigma
              - Partial secondary structure asssigned to this point Slist
              - delta, max distance from the optimal allowed
              - The dynamic programming tables computed W and S
              - A sequence seq of nucleotides {A, U, C, G} of size n
              - A energy model to work with
              - A model of BP considered BPconsidered
       Output: All secondary structure that terminates with a BP in the next element in sigma over seq and at delta from optimal"""
    if sigma == []:
        return [(Slist)]
    (i,j) = sigma.pop()
    
    if j<= i:
        return DeltaBackTrackW(sigma, Slist, delta, W, S, seq, model=model, BPconsidered=BPconsidered)
    
    resu = []
    
    delt= S[i + 1][j - 1] + Energy(seq, i, j, model=model, BPconsidered=BPconsidered) - S[i][j]
    if delta - delt>= 0 and isValid(seq, i, j, BPconsidered=BPconsidered):

        newsigma = sigma + [(i + 1, j - 1)] 

        newSlist = [(i, j)] + Slist
        resu += DeltaBackTrackS(newsigma, newSlist, delta - delt, W, S, seq, model=model, BPconsidered=BPconsidered)
    
    delt= W[i + 1][j - 1] - S[i][j]
    if delta - delt>= 0 and isValid(seq, i, j, BPconsidered=BPconsidered):

        newsigma = sigma + [(i + 1, j - 1)] 
        newSlist = [(i, j)] + Slist
        resu += DeltaBackTrackW(newsigma, newSlist, delta - delt, W, S, seq, model=model, BPconsidered=BPconsidered)
    
    return resu

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

def main_unitary(seq, model="Unitary", BPconsidered="Nussinov", debug = 0, output_format = "png", name_file = "", show=False):
    """Input: - A Sequence seq of letters in {A, U, C, G} of size n
             - A energy model to work with
             - A model of BP considered BPconsidered
             - debug, print for debug
             - output_format is the type of file that we want to save, specify None for no save.
             - name_file is the name of the output Varna file if specified.
             - show is set to False if we want to observe the output directly or not (interactively).
      Output: Build a MFE structure, print it and draw it  in the Unitary model"""
    M = FillMatUnitary(seq, model=model, BPconsidered=BPconsidered)
    if debug:
        print(M)
    S = OptimalBackTrackUnitary(0, len(seq) - 1, M, seq, model=model, BPconsidered=BPconsidered)
    if debug:
        print(S)
    return S#convert(seq, [S], output_format=output_format, name_file=name_file, show=show, seeder=42)


def delta_main_unitary(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = "png", name_file = "", show=False, full = 0):
    #Warning this function is bugged for now
    """Input: - A Sequence seq of letters in {A, U, C, G} of size n
             - A energy model to work with
             - A model of BP considered BPconsidered
             - delta, max distance from the optimal allowed
             - debug, print for debug
             - output_format is the type of file that we want to save, specify None for no save.
             - name_file is the name of the output Varna file if specified.
             - show is set to False if we want to observe the output directly or not (interactively).
      Output: Build a MFE structure, print it and draw it  in the Unitary model"""
    M = FillMatUnitary(seq, model=model, BPconsidered=BPconsidered)
    if debug:
        print(M)
    Slist = DeltaBackTrackUnitary([(0, len(seq) - 1)], [], delta, M, seq, model=model, BPconsidered=BPconsidered)
    for S in Slist:
        if debug:
            print(S)   
    return Slist#convert(seq, Slist, output_format=output_format, name_file=name_file, show=show, full = full)
        

#seq = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAPCUGGAGGUCCUGUGUPCGAUCCACAGAAUUCGCACCA"
#seq = "UAUUUUUAAACGCGCG" #Tried a few variations
#main_unitary(seq, model="Unitary", BPconsidered="Nussinov", debug = 0, output_format = None, name_file = "", show=True)

#delta_main_unitary("AAAUUUUUUAAAAAAUUU", model="Unitary", BPconsidered="Nussinov", delta = 1, debug = 0, output_format = None, name_file = "", show=False)

def main_stacking(seq, model="Unitary", BPconsidered="Nussinov", debug = 0, output_format = "png", name_file = "", show=False):
    """Input: - A Sequence seq of letters in {A, U, C, G} of size n
             - A energy model to work with
             - A model of BP considered BPconsidered
             - debug, print for debug
             - output_format is the type of file that we want to save, specify None for no save.
             - name_file is the name of the output Varna file if specified.
             - show is set to False if we want to observe the output directly or not (interactively).
      Output: Build a MFE structure, print it and draw it  in the Unitary model"""
    W, S = FillMatStacking(seq, BPconsidered=BPconsidered)
    if debug:
        print(W)
        print(S)
    St = OptimalBackTrackW(0, len(seq) - 1, W, S, seq, model=model, BPconsidered=BPconsidered)
    if debug:
        print(St)
    convert(seq, St, output_format=output_format, name_file=name_file, show=show)

def main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", debug = 0, output_format = "png", name_file = "", show=False):
    """Input: - A Sequence seq of letters in {A, U, C, G} of size n
             - A energy model to work with
             - A model of BP considered BPconsidered
             - debug, print for debug
             - output_format is the type of file that we want to save, specify None for no save.
             - name_file is the name of the output Varna file if specified.
             - show is set to False if we want to observe the output directly or not (interactively).
      Output: Build a MFE structure, print it and draw it  in the Unitary model"""
    E, S = FillMatStacking2(seq, BPconsidered=BPconsidered)
    if debug:
        print(E)
        print(S)
    St = OptimalBackTrackE2(0, len(seq) - 1, E, S, seq, model=model, BPconsidered=BPconsidered)
    if debug:
        print(St)
    convert(seq, [St], output_format=output_format, name_file=name_file, show=show, seeder=42)
    
def delta_main_stacking(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = "png", name_file = "", show=False):
    """Input: - A Sequence seq of letters in {A, U, C, G} of size n
             - A energy model to work with
             - A model of BP considered BPconsidered
             - delta, max distance from the optimal allowed
             - debug, print for debug
             - output_format is the type of file that we want to save, specify None for no save.
             - name_file is the name of the output Varna file if specified.
             - show is set to False if we want to observe the output directly or not (interactively).
      Output: Build a MFE structure, print it and draw it  in the Unitary model"""
    W, S = FillMatStacking(seq, BPconsidered=BPconsidered)
    if debug:
        print(W)
        print(S)
    Slist = DeltaBackTrackW([(0, len(seq) - 1)], [], delta, W, S, seq, model=model, BPconsidered=BPconsidered)
    for St in Slist:
        if debug:
            print(St)
    convert(seq, Slist, output_format=output_format, name_file=name_file, show=show)

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
    return Slist#convert(seq, Slist, output_format=output_format, name_file=name_file, show=show)
         
#"AAAUUUUUUAAGCAAUUU"
#main_stacking("AAAUUUUUUAAGCAAUUU", model="Unitary", BPconsidered="Nussinov", debug = 1, output_format = None, name_file = "", show=False)


#main_stacking("AAGGGAUUUACCCGGGUU", model="Unitary", BPconsidered="Nussinov", debug = 1, output_format = None, name_file = "", show=False)

#main_stacking("AAGGGAUUUACCCGGGUU", model="Unitary", BPconsidered="Nussinov", debug = 1, output_format = None, name_file = "", show=False)
#delta_main_stacking("AAUUUUGG", model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 1, output_format = None, name_file = "", show=False)
#main_stacking("AAGGGAUUUACCCGGGUU", model="Unitary", BPconsidered="Nussinov", debug = 1, output_format = None, name_file = "", show=False)
#print("START")
#delta_main_stacking("AAGGGAUUUACCCGGGUU", model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 1, output_format = None, name_file = "", show=False)
#print("STOP")
#main_stacking("AAGGGAUUUACCCGGGUU", model="Unitary", BPconsidered="Nussinov", debug = 1, output_format = None, name_file = "", show=False)

#print("START")
#delta_main_stacking("AAGGGAUUUACCCGGGUUAAGGGAUUUACCCGGGUU", model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 1, output_format = None, name_file = "", show=False)
#print("STOP")
#main_stacking2("AAGGGAUUUACCCGGGUU", model="Unitary", BPconsidered="Nussinov", debug = 1, output_format = None, name_file = "", show=False)

#print("START")
#delta_main_stacking2("AAGGGAUUUACCCGGGUU", model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 1, output_format = None, name_file = "", show=False)
#print("STOP")

#main_stacking2("AAGGGAUUUACCCGGGUU", model="Unitary", BPconsidered="Nussinov", debug = 1, output_format = None, name_file = "", show=False)
#delta_main_stacking2("GAGCUCCCCGGGGUGCAC", model="Unitary", BPconsidered="Nussinov", delta = 1, debug = 1, output_format = None, name_file = "", show=True)

"""
#h1 = 3, h2 = 3
for h2 in range(2, 15):
    seq = "GAGCUC" + "C"*h2 + "G"*h2 + "GUGCAC"
    nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
    if nbseq != 1:
        print("ALERT", seq)
        
#h1 >= 4, h2 = 3
for h2 in range(2, 15):
    seq = "GAGGCCUC" + "C"*h2 + "G"*h2 + "GUGCAC"
    nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
    if nbseq != 1:
        print("ALERT", seq)

#h1 >= 4 ou 5, h2 = 3
for h1 in range(4, 15):
    for h2 in range(2, 15):
        seq = "GA" + "G"*(h1 - 2) + "C"*(h1 - 2) + "UC" + "C"*h2 + "G"*h2 + "GUGCAC"
        nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
        if nbseq != 1:
            print("ALERGAGGCCUC" + "C"*h2 + "G"*h2 + "GUGCAC"
    nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
    if nbseq != 1:
        print("ALERT", seq)

#h1 >= 4 ou 5, h2 = 3
for h1 in range(4, 15):
    for h2 in range(2, 15):
        seq = "GA" + "G"*(h1 - 2) + "C"*(h1 - 2) + "UC" + "C"*h2 + "G"*h2 + "GUGCAC"
        nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
        if nbseq != 1:
            print("ALERT", seq)
        
#h1 >= 4 ou 5, h2 = 4
for h1 in range(4, 15):
    for h2 in range(2, 15):
        seq = "GA" + "G"*(h1 - 2) + "C"*(h1 - 2) + "UC" + "C"*h2 + "G"*h2 + "GUGGCCAC"
        
        nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
        if nbseq != 1:
            print("ALERT1", seq)

#Issues with "GAGGGCCCUCCCGGGUCCGGAC" "GAGGCCUCCCGGGUCCGGAC","GAGGGGCCCCUCCCGGGUCCGGAC"
        
#h1 >= 4 ou 5, h2 >= 5
for h1 in range(4, 15):
    for h2 in range(2, 15):
        for h3 in range(5, 15):
            seq = "GA" + "G"*(h1 - 2) + "C"*(h1 - 2) + "UC" + "C"*h2 + "G"*h2 + "GUG" + "C"*(h3 - 3) + "G"*(h3 - 3) + "CAC"
        
        nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
        if nbseq != 1:
            print("ALERT2", seq):
    for h2 in range(2, 15):
        seq = "GA" + "G"*(h1 - 2) + "C"*(h1 - 2) + "UC" + "C"*h2 + "G"*h2 + "GUGGCCAC"
        
        nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
        if nbseq != 1:
            print("ALERT1", seq)

#Issues with "GAGGGCCCUCCCGGGUCCGGAC" "GAGGCCUCCCGGGUCCGGAC","GAGGGGCCCCUCCCGGGUCCGGAC"
        
#h1 >= 4 ou 5, h2 >= 5
for h1 in range(4, 15):
    for h2 in range(2, 15):
        for h3 in range(5, 15):
            seq = "GA" + "G"*(h1 - 2) + "C"*(h1 - 2) + "UC" + "C"*h2 + "G"*h2 + "GUG" + "C"*(h3 - 3) + "G"*(h3 - 3) + "CAC"
        
        nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
        if nbseq != 1:
            print("ALERT2", seq)

#h1 >= 4 ou 5, h2 >= 5
for h1 in range(4, 11):
    for h2 in range(2, 15):
        for h3 in range(h1, 11):
            seq = "GGGA" + "G"*(h1 - 4) + "C"*(h1 - 4) + "U" + "C"*(h2 + 3) + "G"*(h2 + 3) + "U" + "G"*(h3 - 4) + "C"*(h3 - 4) + "ACCC"
            nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
            if nbseq != 1:
                print("ALERT3", h1, h2, h3, seq)
"""
              

#for h1 in range(4, 11):
#    for h2 in range(2, 15):
#        for h3 in range(h1, 11):
#            seq = "GGA" + "G"*(h1 - 3) + "C"*(h1 - 3) + "U" + "C"*(h2 + 2) + "G"*(h2 + 2) + "U" + "G"*(h3 - 3) + "C"*(h3 - 3) + "ACC"
#            nbseq = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
#            if nbseq != 1:
#                print("ALERT", h1, h2, h3, seq)
