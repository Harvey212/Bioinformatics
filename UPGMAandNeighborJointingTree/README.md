# hw5

<葉冠宏 + 108753208>

## Description

* Write R or Python script to build a unweighted pair group method with arithmetic mean (UPGMA) or Neighbor Jointing tree given a distance matrix.
* Creating your own script with your student ID, ie. hw5_105753026.R or hw5_105753026.py.

## Files

* distanceMatrix1.txt where the label of leaves start from A to Z, i.e., A, B, C, D for four taxa
* UPGMA-tree.nwk, UPGMA tree in newick format with branch length
* NJ-tree.nwk, NJ tree in newick format with branch length

## Parameters

* input: a given distance matrix file (ex. matrix.txt)
* method: UPGMA | NJ
* output: the output tree file

## Command

Executing your code with the following command.

```R
Rscript hw5_studentID.R --input distanceMatrix1.txt --method UPGMA --output UPGMA-tree.nwk
Rscript hw5_studentID.R --input distanceMatrix1.txt --method NJ --output NJ-tree.nwk
```

```Python
python3 hw5_studentID.py --input distanceMatrix1.txt --method UPGMA --output UPGMA-tree.nwk
python3 hw5_studentID.py --input distanceMatrix1.txt --method NJ --output NJ-tree.nwk
```

## Evaluation

10 testing data

Correct answer gets 10 points of each testing data.

### Penalty

* High code similarity to others: YOUR SCORE = 0
