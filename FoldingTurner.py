import os 
def fold_turner(seq):
    """
    Input:
        * seq, a sequence that we designed
        * target_struct, the secondary structure that we target 
        * h, a maximal BP distance 
        * pr, a boolean if we want some info printed or not
    Output: 
        * The corresponding structure in Turner
    """
    with open("designedseq.fa", "w") as f:
        f.write(">sequence1\n")
        f.write(seq + "\n")
    os.system("RNAsubopt -v -s -d2 -e 0 < designedseq.fa > subopts.txt")
    with open("subopts.txt", "r") as f:
        li = f.readlines()[2:]  
    os.system("RNAfold -p -d2 --noLP < designedseq.fa > designedseqandturner.out")

    Turner_struct = (li[0].strip().split(" "))[0]
    #print(li)
    nb = len(li)
    os.system("rm designedseq.fa")
    os.system("rm subopts.txt")          
    return Turner_struct, nb
