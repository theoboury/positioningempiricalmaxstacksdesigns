# positioningempiricalmaxstacksdesigns
# Copyright (C) 2026 THEO BOURY 

from checkSeparability import fullSeparable, dbn_to_tree, filter, ssparse, first_modulo_separable
from FoldingTurner import fold_turner
from foldingBPandStacking import main_unitary, delta_main_stacking2 #Warning delta_main_unitary is bugged for now
from foldingBP import main_unitary_only_one
from SecondaryStructureGeneration import ssrandom #ssrandom(n,count,count_stacked,theta,min_helix)
from RandomCompatible import choose_random_seq
import csv
import random



def sstopairs(ss):
    """
    Input:
        * ss, a secondary structure as a list of complementary base pairs
    Output: 
        * The list of all base pairs
    """
    if ss is None:
        return None
    res = []
    for i in range(len(ss)):
        if ss[i]>i:
            j = ss[i]
            res.append((i,j))
    return set(res)

def pairstoss(n, li):
    resu = ["." for i in range(n)]
    for (i,j) in li:
        resu[i] = "("
        resu[j] = ")"

    return "".join(resu)

def create_stats_from_Stacking_A_only_nom3o_nom5(n, iteration,theta,min_helix, restart=1,last_index=-1):
    c = 'a'
    if restart:
        c = 'w'
    with open('ResultsfromStacking.csv', c, newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        if restart:
            elem = ["ss", "seq", "Separable(Aonly)", "BPDesign", "BPfold", "TurnerDesign", "Turnerfold"]
            writer.writerow(elem)
        count,count_stacked = {}, {}
        for i in range(last_index+1, iteration):
            resu = []
            ss = ssrandom(n,count,count_stacked,theta,min_helix)
            t = dbn_to_tree(ssparse(ss))
            while not filter(t):
                ss = ssrandom(n,count,count_stacked,theta,min_helix)
                t = dbn_to_tree(ssparse(ss))
            print("iteration", i, " ss", ss)
            seq = choose_random_seq(t, withA=True)
            timeout = 0
            Slist  = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
            Slist0struct = pairstoss(n, Slist[0])
            while ((len(Slist) != 1) or  (Slist0struct != ss)):
                seq = choose_random_seq(t, withA=True)
                Slist  = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
                
                if timeout >= 1000:
                    ss = ssrandom(n,count,count_stacked,theta,min_helix)
                    t = dbn_to_tree(ssparse(ss))
                    while not filter(t):
                        ss = ssrandom(n,count,count_stacked,theta,min_helix)
                        t = dbn_to_tree(ssparse(ss))
                    print("iteration", i, " again, ss", ss)
                    seq = choose_random_seq(t, withA=True)
                    timeout = 0
                else:
                    timeout+=1
                Slist0struct = pairstoss(n, Slist[0])
            resu.append(ss)
            resu.append(seq)
            #print("seq", seq, "\n")
            Separable = "False"
            if fullSeparable(seq, ss):
                Separable = "True"
            resu.append(Separable)
            BPstruct0, nb  = main_unitary_only_one(seq, model="Unitary", BPconsidered="Nussinov")
            #print(nb)
            BPDesign= "False"
            BPstruct = pairstoss(n, BPstruct0)
            if nb == 1 and BPstruct == ss:
                BPDesign = "True"
            resu.append(BPDesign)
            resu.append(BPstruct)
            Turnerss, nb2 = fold_turner(seq)
            TurnerDesign = "False"
            if nb2 == 1 and Turnerss == ss:
                TurnerDesign = "True"
            resu.append(TurnerDesign)
            resu.append(Turnerss)

            writer.writerow(resu)


def from_stacking_read_stats_from_csv(name):
    isTurnernotBP = 0
    isTurnerandBP = 0
    isBPnotTurner = 0
    isnotBPnotTurner = 0

    isTurnernotSeparable=0
    isTurnerandSeparable=0
    isSeparablenotTurner=0
    isnotSeparablenotTurner=0
    with open(name, 'r') as csvfile:
        for line in csvfile.readlines()[1:]:
            _, _, Separable, BPDesign, _, TurnerDesign, _ = line.split(' ')
            if TurnerDesign == "True" and BPDesign == "False":
                isTurnernotBP+=1
            if TurnerDesign == "True" and BPDesign == "True":
                isTurnerandBP+=1
            if TurnerDesign == "False" and BPDesign == "True":
                isBPnotTurner+=1
            if TurnerDesign == "False" and BPDesign == "False":
                isnotBPnotTurner+=1

            if TurnerDesign == "True" and Separable == "False":
                isTurnernotSeparable+=1
            if TurnerDesign == "True" and Separable == "True":
                isTurnerandSeparable+=1
            if TurnerDesign == "False" and Separable == "True":
                isSeparablenotTurner+=1
            if TurnerDesign == "False" and Separable == "False":
                isnotSeparablenotTurner+=1
    print("isTurnernotBP:", isTurnernotBP," isTurnerandBP:", isTurnerandBP, " isBPnotTurner:", isBPnotTurner, "isnotBPnotTurner:", isnotBPnotTurner, "\n")
    print("isTurnernotSeparable:", isTurnernotSeparable," isTurnerandSeparable:", isTurnerandSeparable, " isSeparablenotTurner:", isSeparablenotTurner, "isnotSeparablenotTurner:", isnotSeparablenotTurner, "\n")


def create_stats_from_Separable_A_only_nom3o_nom5(n, iteration,theta,min_helix, restart=1,last_index=-1):
    c = 'a'
    if restart:
        c = 'w'
    with open('ResultsfromSeparable.csv', c, newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        if restart:
            elem = ["ss", "seq", "StackingDesign","StackingFold", "TurnerDesign", "Turnerfold"]
            writer.writerow(elem)
        count,count_stacked = {}, {}



        for i in range(last_index+1, iteration):
            resu = []
            ss = ssrandom(n,count,count_stacked,theta,min_helix)
            t = dbn_to_tree(ssparse(ss))
            while not filter(t):
                ss = ssrandom(n,count,count_stacked,theta,min_helix)
                t = dbn_to_tree(ssparse(ss))
            print("iteration", i, " ss", ss)
            #ss_useful = sstopairs(ssparse(ss))
            #print("ss", ss_useful)
            seq = None
            while seq is None:
                seq = first_modulo_separable(t, modulolimit=2)
            #seq = choose_random_seq(t, withA=True)
            resu.append(ss)
            resu.append(seq)
            Slist  = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
            StackDesign="False"
            Stackstruct = pairstoss(n, Slist[0])
            if len(Slist) == 1 and Stackstruct == ss:
                StackDesign="True"
            resu.append(StackDesign)
            resu.append(Stackstruct)
            Turnerss, nb2 = fold_turner(seq)
            TurnerDesign = "False"
            if nb2 == 1 and Turnerss == ss:
                TurnerDesign = "True"
            resu.append(TurnerDesign)
            resu.append(Turnerss)
            writer.writerow(resu)



def from_separable_read_stats_from_csv(name):
    isTurnernotStacking = 0
    isTurnerandStacking = 0
    isStackingnotTurner = 0
    isnotStackingnotTurner = 0
    with open(name, 'r') as csvfile:
        for line in csvfile.readlines()[1:]:
            _, _, StackDesign, _, TurnerDesign, _ = line.split(' ')
            if TurnerDesign == "True" and StackDesign == "False":
                isTurnernotStacking+=1
            if TurnerDesign == "True" and StackDesign == "True":
                isTurnerandStacking+=1
            if TurnerDesign == "False" and StackDesign == "True":
                isStackingnotTurner+=1
            if TurnerDesign == "False" and StackDesign == "False":
                isnotStackingnotTurner+=1
    print("isTurnernotStacking:", isTurnernotStacking," isTurnerandStacking:", isTurnerandStacking, " isStackingnotTurner:", isStackingnotTurner, "isnotStackingnotTurner:", isnotStackingnotTurner, "\n")

def create_stats_from_Stacking_A_only_withm3o_withm5(n, iteration,theta,min_helix, restart=1,last_index=-1):
    c = 'a'
    if restart:
        c = 'w'
    with open('ResultsfromStackingwithm3oandm5.csv', c, newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        if restart:
            elem = ["ss", "seq", "TurnerDesign", "Turnerfold"]
            writer.writerow(elem)
        count,count_stacked = {}, {}
        for i in range(last_index+1, iteration):
            resu = []
            ss = ssrandom(n,count,count_stacked,theta,min_helix)
            t = dbn_to_tree(ssparse(ss))
            while filter(t):
                ss = ssrandom(n,count,count_stacked,theta,min_helix)
                t = dbn_to_tree(ssparse(ss))
            print("iteration", i, " ss", ss)
            seq = choose_random_seq(t, withA=True)
            timeout = 0
            Slist  = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
            Slist0struct = pairstoss(n, Slist[0])
            while ((len(Slist) != 1) or  (Slist0struct != ss)):
                seq = choose_random_seq(t, withA=True)
                Slist  = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
                
                if timeout >= 1000:
                    ss = ssrandom(n,count,count_stacked,theta,min_helix)
                    t = dbn_to_tree(ssparse(ss))
                    while filter(t):
                        ss = ssrandom(n,count,count_stacked,theta,min_helix)
                        t = dbn_to_tree(ssparse(ss))
                    print("iteration", i, " again, ss", ss)
                    seq = choose_random_seq(t, withA=True)
                    timeout = 0
                else:
                    timeout+=1
                Slist0struct = pairstoss(n, Slist[0])
            resu.append(ss)
            resu.append(seq)
            Turnerss, nb2 = fold_turner(seq)
            TurnerDesign = "False"
            if nb2 == 1 and Turnerss == ss:
                TurnerDesign = "True"
            resu.append(TurnerDesign)
            resu.append(Turnerss)

            writer.writerow(resu)


def from_stacking_withm3om5_read_stats_from_csv(name):
    isTurner = 0
    with open(name, 'r') as csvfile:
        for line in csvfile.readlines()[1:]:
            _, _, TurnerDesign, _ = line.split(' ')
            if TurnerDesign == "True":
                isTurner+=1
    print("isTurner:", isTurner, "\n")


def refine_stats_from_Stacking_A_only_withm3o_withm5():
    with open('ResultsfromStackingwithm3oandm5increased.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        elem = ["ss", "seq", "TurnerDesign", "Turnerfold", "random_compatible_seq", "random_compatible_TurnerDesign", "random_compatible_Turnerfold"]
        writer.writerow(elem)
        i = -1
        with open('ResultsfromStackingwithm3oandm5.csv', 'r') as readfile:
            for line in readfile.readlines()[1:]:
                i+=1
                ss, seq, TurnerDesign, Turnerfold = line.strip().split(' ')
                print("iteration", i, " ss", ss)
                resu = [ss, seq, TurnerDesign, Turnerfold]
                t = dbn_to_tree(ssparse(ss))
                random_compatible_seq = choose_random_seq(t, withA=True)
                random_compatible_Turnerfold, nb2 = fold_turner(random_compatible_seq)
                random_compatible_TurnerDesign = "False"
                if nb2 == 1 and random_compatible_Turnerfold == ss:
                    random_compatible_TurnerDesign = "True"
                resu.append(random_compatible_seq)
                resu.append(random_compatible_TurnerDesign)
                resu.append(random_compatible_Turnerfold)
                print(resu, "\n")
                writer.writerow(resu)

def from_stacking_withm3om5increased_read_stats_from_csv(name):
    isTurnerandrandom_compatible_TurnerDesign = 0
    isTurnernotrandom_compatible_TurnerDesign = 0
    israndom_compatible_TurnerDesignnotTurner = 0
    isnotrandom_compatible_TurnerDesignnotTurner = 0
    with open(name, 'r') as csvfile:
        for line in csvfile.readlines()[1:]:
            _, _, TurnerDesign, _, _, random_compatible_TurnerDesign, _ = line.split(' ')
            if TurnerDesign == "True" and random_compatible_TurnerDesign == "True":
                isTurnerandrandom_compatible_TurnerDesign+=1
            if TurnerDesign == "True" and random_compatible_TurnerDesign == "False":
                isTurnernotrandom_compatible_TurnerDesign+=1
            if TurnerDesign == "False" and random_compatible_TurnerDesign == "True":
                israndom_compatible_TurnerDesignnotTurner+=1
            if TurnerDesign == "False" and random_compatible_TurnerDesign == "False":
                isnotrandom_compatible_TurnerDesignnotTurner+=1
    print("isTurnerandrandom_compatible_TurnerDesign:", isTurnerandrandom_compatible_TurnerDesign, "\n",
          "isTurnernotrandom_compatible_TurnerDesign",isTurnernotrandom_compatible_TurnerDesign, "\n",
          "israndom_compatible_TurnerDesignnotTurner", israndom_compatible_TurnerDesignnotTurner, "\n",
          "isnotrandom_compatible_TurnerDesignnotTurner", isnotrandom_compatible_TurnerDesignnotTurner
          , "\n")



def stacking_vs_BP_A_only_nom3o_nom5(n, iteration,theta,min_helix, restart=1,last_index=-1):
    c = 'a'
    if restart:
        c = 'w'
    with open('ResultsStackingvsBP.csv', c, newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        if restart:
            elem = ["ss", "Stacking_seq", "Stacking_TurnerFold", "Stacking_TurnerDesign", "BP_seq", "BP_TurnerFold", "BP_TurnerDesign", "nb_it_more_for_finding_BP"]
            writer.writerow(elem)
        count,count_stacked = {}, {}
        for i in range(last_index+1, iteration):
            resu = []
            ss = ssrandom(n,count,count_stacked,theta,min_helix)
            t = dbn_to_tree(ssparse(ss))
            while not filter(t):
                ss = ssrandom(n,count,count_stacked,theta,min_helix)
                t = dbn_to_tree(ssparse(ss))
            print("iteration", i, " ss", ss)
            timeout = 0
            BPDesign= False
            StackingDesign= False
            while (not BPDesign) or (not StackingDesign):
                #(len(Slist) != 1) or  (Slist0struct != ss)
                seq = choose_random_seq(t, withA=True)
                if not StackingDesign:
                    Slist  = delta_main_stacking2(seq, model="Unitary", BPconsidered="Nussinov", delta = 0, debug = 0, output_format = None, name_file = "", show=False)
                    Slist0struct = pairstoss(n, Slist[0])
                    if (len(Slist) == 1) and (Slist0struct == ss):
                        StackingDesign = True
                        Stacking_seq = seq
                        it_stack=timeout
                        #print("it_stack", it_stack, "\n")
                if not BPDesign:
                    #BPstruct  = main_unitary(seq, model="Unitary", BPconsidered="Nussinov", debug = 0, output_format = None, name_file = "", show=False)
                    BPstruct, nb  = main_unitary_only_one(seq, model="Unitary", BPconsidered="Nussinov")
                    BPstruct = pairstoss(n, BPstruct)
                    if nb == 1 and BPstruct == ss:
                        BPDesign = True
                        BP_seq = seq
                        it_BP=timeout
                        #print("it_BP", it_BP, "\n")
                if timeout >= 10000:
                    print("nb", nb)
                    print("BPDesign", BPDesign,  "StackingDesign", StackingDesign)
                    ss = ssrandom(n,count,count_stacked,theta,min_helix)
                    t = dbn_to_tree(ssparse(ss))
                    while not filter(t):
                        ss = ssrandom(n,count,count_stacked,theta,min_helix)
                        t = dbn_to_tree(ssparse(ss))
                    print("iteration", i, " again, ss", ss)
                    timeout = 0
                else:
                    timeout+=1



            Stacking_TurnerFold, nbStack = fold_turner(Stacking_seq)
            BP_TurnerFold, nbBP = fold_turner(BP_seq)
            Stacking_TurnerDesign = "False"
            BP_TurnerDesign = "False"
            if nbStack == 1 and Stacking_TurnerFold == ss:
                Stacking_TurnerDesign = "True"
            if nbBP == 1 and BP_TurnerFold == ss:
                BP_TurnerDesign = "True"
            nb_it_more_for_finding_BP = it_BP - it_stack
            resu = [ss, Stacking_seq, Stacking_TurnerFold, Stacking_TurnerDesign, BP_seq, BP_TurnerFold, BP_TurnerDesign, nb_it_more_for_finding_BP]
            writer.writerow(resu)


def stacking_vs_BP_read_stats_from_csv(name):
    isStacking_TurnerDesignandBP_TurnerDesign = 0
    isStacking_TurnerDesignnotBP_TurnerDesign = 0
    isBP_TurnerDesignnotStacking_TurnerDesign = 0
    isnotBP_TurnerDesignnotStacking_TurnerDesign = 0


    with open(name, 'r') as csvfile:
        res_it = 0
        num = 0
        for line in csvfile.readlines()[1:]:
            num+=1
            _, _, _, Stacking_TurnerDesign, _, _, BP_TurnerDesign, nb_it_more_for_finding_BP = line.strip().split(' ')
            res_it+= float(nb_it_more_for_finding_BP)
            if Stacking_TurnerDesign == "True" and BP_TurnerDesign == "True":
                isStacking_TurnerDesignandBP_TurnerDesign+=1
            if Stacking_TurnerDesign == "True" and BP_TurnerDesign == "False":
                isStacking_TurnerDesignnotBP_TurnerDesign+=1
            if Stacking_TurnerDesign == "False" and BP_TurnerDesign == "True":
                isBP_TurnerDesignnotStacking_TurnerDesign+=1
            if Stacking_TurnerDesign == "False" and BP_TurnerDesign == "False":
                isnotBP_TurnerDesignnotStacking_TurnerDesign+=1
    print("isStacking_TurnerDesignandBP_TurnerDesign:", isStacking_TurnerDesignandBP_TurnerDesign, "\n",
    "isStacking_TurnerDesignnotBP_TurnerDesign:", isStacking_TurnerDesignnotBP_TurnerDesign, "\n",
    "isBP_TurnerDesignnotStacking_TurnerDesign:", isBP_TurnerDesignnotStacking_TurnerDesign, "\n",
    "isnotBP_TurnerDesignnotStacking_TurnerDesign:", isnotBP_TurnerDesignnotStacking_TurnerDesign, "\n",
    "In mean, nb of additional iteration necessary to get a design in BP (not taking restart into account)", res_it/num)

"""Examples of code usages:
To begin:
create_stats_from_Stacking_A_only_nom3o_nom5(n=150, iteration=2000,theta=3,min_helix=3, restart=1,last_index=-1)
To restart:
create_stats_from_Stacking_A_only_nom3o_nom5(n=150, iteration=2000,theta=3,min_helix=3, restart=0,last_index=lastcomputed)
To see proportions:
from_stacking_read_stats_from_csv('ResultsfromStacking.csv')

To begin:
create_stats_from_Separable_A_only_nom3o_nom5(n=150, iteration=2000,theta=3,min_helix=3, restart=1,last_index=-1)
To restart:
create_stats_from_Separable_A_only_nom3o_nom5(n=150, iteration=2000,theta=3,min_helix=3, restart=0,last_index=lastcomputed)
To see proportions:
from_separable_read_stats_from_csv('ResultsfromSeparable.csv')

To begin:
create_stats_from_Stacking_A_only_withm3o_withm5(n=150, iteration=2000,theta=3,min_helix=3, restart=1,last_index=-1)
To restart:
create_stats_from_Stacking_A_only_withm3o_withm5(n=150, iteration=2000,theta=3,min_helix=3, restart=0,last_index=lastcomputed)
To see proportions:
from_stacking_withm3om5_read_stats_from_csv('ResultsfromStackingwithm3oandm5.csv')

To begin:
refine_stats_from_Stacking_A_only_withm3o_withm5()
To see proportions:
from_stacking_withm3om5increased_read_stats_from_csv('ResultsfromStackingwithm3oandm5increased.csv')

To begin:
stacking_vs_BP_A_only_nom3o_nom5(n=50, iteration=2000,theta=3,min_helix=3, restart=1,last_index=-1)
To restart:
stacking_vs_BP_A_only_nom3o_nom5(n=50, iteration=2000,theta=3,min_helix=3, restart=0,last_index=lastcomputed)
To see proportions:
stacking_vs_BP_read_stats_from_csv('ResultsStackingvsBP.csv')
"""



