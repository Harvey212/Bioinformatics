# hw2
<葉冠宏 + 108753208>
## Description

* Write R or Python script to calculate the sum-of-pair score (SoP) of the multiple sequence alignment.
* Creating your own script with your student ID, ie. hw2_105753026.R.
* In this program, library Biostrings is only used to parse fasta file.

## File

* hw2_ref.R: You can start from this reference code, and try to write your own comment in English
* pam100.txt
* pam250.txt
* test.fasta

## Parameters

* input: fasta file (ex. test.fasta)
* score: score file (ex. pam250.txt)
* gopen: gap open penalty
* gextend: gap extend penalty

## Command

Executing your code with the following command.

```R
Rscript hw2_studentID.R --input test.fasta --score pam250.txt --gopen -10 --gextend -2
```

```Python
python3 hw2_studentID.R --input test.fasta --score pam250.txt --gopen -10 --gextend -2
```
The answer is 1047. You should print it on the screen.

## Evaluation

10 testing data

```R
Rscript hw2_studentID.R --input test.fasta --score pam250.txt --gopen -10 --gextend -2
Rscript hw2_studentID.R --input test2.fasta --score pam100.txt --gopen -10 --gextend -2
Rscript hw2_studentID.R --input data/test3.fasta --score pam/pam1.txt --gopen -10 --gextend -2
```

```Python
python3 hw2_studentID.R --input test.fasta --score pam250.txt --gopen -10 --gextend -2
python3 hw2_studentID.R --input test2.fasta --score pam100.txt --gopen -10 --gextend -2
python3 hw2_studentID.R --input data/test3.fasta --score pam/pam1.txt --gopen -10 --gextend -2
```


Correct answer gets 10 points of each testing data.

### Penalty

* High code similarity to others: YOUR SCORE = 0

