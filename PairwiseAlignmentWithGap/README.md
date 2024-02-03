# hw4
<葉冠宏 108753208>
## Description

* Write R or Python script to perform pairwise alignment.
* You should divide gap scores into two types(open gap, extend gap) based on your hw3.
* Creating your own script with your student ID, ie. hw4_105753026.R or hw4_105753026.py.
* In this program, library Biostrings is only used to parse fasta file.

## Files

* hw4_ref.R: You can start from this reference code, and try to write your own comment in English.
* read_save_fasta.R: An example of reading and saving fasta file with Biostrings library.
* pam100.txt
* pam250.txt
* test.fasta
* result.fasta: An example of output file in FASTA format.

## Parameters

* input: fasta file (ex. test.fasta)
* score: score file (ex. pam250.txt)
* aln: global|local
* gap_open: open gap score
* gap_extend: extend gap score
* output: .fasta file

## Command

Executing your code with the following command.

```R
Rscript hw4_studentID.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta
```

```Python
python3 hw4_studentID.R --input test1.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result1.fasta
```
You should output your answer into fasta format.

## Evaluation

10 testing data

```R
Rscript hw4_studentID.R --input test1.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result1.fasta
Rscript hw4_studentID.R --input test2.fasta --score pam100.txt --aln local --gap_open -8 --gap_extend -5 --output result2.fasta
Rscript hw4_studentID.R --input data/test3.fasta --score pam/pam100.txt --aln local --gap_open -5 --gap_extend -3 --output output/result3.fasta
```

```Python
python3 hw4_studentID.R --input test1.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result1.fasta
python3 hw4_studentID.R --input test2.fasta --score pam100.txt --aln local --gap_open -8 --gap_extend -5 --output result2.fasta
python3 hw4_studentID.R --input data/test3.fasta --score pam/pam100.txt --aln local --gap_open -5 --gap_extend -3 --output output/result3.fasta
```

Correct answer gets 10 points of each testing data.

### Penalty

* High code similarity to others: YOUR SCORE = 0

