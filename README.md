# positioningempiricalmaxstacksdesigns

Copyright (C) 2026 THEO BOURY

This git contains the functions to reproduce the experiments to position the maxStacks energy model for RNA Inverse Folding with respect to other energy models as depicted in #LINKTOHALTOINSERTWHENAVAILABLE. All codes are in Python3.

### Dependencies

- Mandatory Python libraries:
    - random
    - math
    - csv
    - argparse
- Outside library:
	- ViennaRNA-2.5.1 (see https://www.tbi.univie.ac.at/RNA/ for more details.)


#### Files and repositories

The different files are:
    checkSeparability.py: check the separability criterion
    foldingBP/Stacking/Turner.py: compute the fold  structure from a sequence in the maxBP/maxStacks/Turner energy model
    RandomCompatible.py: create random compatible sequences
    SecondaryStructureGeneration.py: create random structures
	positioningempiricalmaxstacksdesigns.py: a parser to launch the experiments in the command line.

### Launch

All experiments can be launched with the following command: 
```bash
python3 MaxStacksPositioning -n 150 -o FromStackingNoLargeLoop
```
With
- n, the size of sequences to sample
- o is the desired experiments to launch that can be:
    * FromStackingNoLargeLoop: Structures with no large loops. Sequence are randomly sampled stacking designs with A at the unpaired positions. The results are put in 'ResultsfromStacking.csv'.
    * FromSeparableNoLargeLoop: Structures with no large loops. Sequence are randomly sampled separable designs with A at the unpaired positions. The results are put in 'ResultsfromSeparable.csv'.
    * FromStackingOnlyLargeLoop: Structures necessarily with large loops. Sequence are randomly sampled stacking designs with A at the unpaired positions. The results are put in 'ResultsfromStackingwithm3oandm5.csv' and 'ResultsfromStackingwithm3oandm5increased.csv'.
    * StackingVsBP: Structures necessarily with no large loops. Sequence are randomly sampled maxStacks and maxBP designs with A at the unpaired positions. The results are put in 'ResultsStackingvsBP.csv'.

### Contributors

    Théo Boury




