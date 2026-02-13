# positioningempiricalmaxstacksdesigns
# Copyright (C) 2026 THEO BOURY 

import argparse
from createrandomsequencesandfold import *

parser = argparse.ArgumentParser(prog='MaxStacksPositioning')


parser.add_argument('-n', '--n', type=int, required=True, help="The size of sequences to sample")
parser.add_argument('-e', '--experiment', type=str, required=False, help="The type of experiment to launch, default is FromStackingNoLargeLoop")
args = parser.parse_args()
e = 0
if args.experiment == "FromStackingNoLargeLoop":
    e = 0
elif args.experiment == "FromSeparableNoLargeLoop":
    e = 1
elif args.experiment == "FromStackingOnlyLargeLoop":
    e = 2
elif args.experiment == "StackingVsBP":
    e = 3

#structure = "((((((((((((((((((((....))))))))((((((((((((....))))))))((((((((((....))))))))))((((((.........))))))))))(((((.........)))))))))))))))))(((((((((((((((((((.........)))))((((((((.........))))))))((((((((....))))))))))))((((((((....))))))))((((((((((....))))))))))))))))))))"

if e == 0:
    create_stats_from_Stacking_A_only_nom3o_nom5(n=args.n, iteration=2000,theta=3,min_helix=3, restart=1,last_index=-1)
    from_stacking_read_stats_from_csv('ResultsfromStacking.csv')
elif e == 1:
    create_stats_from_Separable_A_only_nom3o_nom5(n=args.n, iteration=2000,theta=3,min_helix=3, restart=1,last_index=-1)
    from_separable_read_stats_from_csv('ResultsfromSeparable.csv')
elif e == 2:
    create_stats_from_Stacking_A_only_withm3o_withm5(n=args.n, iteration=2000,theta=3,min_helix=3, restart=1,last_index=-1)
    from_stacking_withm3om5_read_stats_from_csv('ResultsfromStackingwithm3oandm5.csv')
    refine_stats_from_Stacking_A_only_withm3o_withm5()
    from_stacking_withm3om5increased_read_stats_from_csv('ResultsfromStackingwithm3oandm5increased.csv')
elif e == 3:
    stacking_vs_BP_A_only_nom3o_nom5(n=50, iteration=args.n,theta=3,min_helix=3, restart=1,last_index=-1)
    stacking_vs_BP_read_stats_from_csv('ResultsStackingvsBP.csv')




