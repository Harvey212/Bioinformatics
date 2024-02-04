# hw3. predict protein subcellular localization

![PredictProtein](/images/img1.png)

### Name: [葉冠宏]
### Student ID: [108753208]

## Description
Perform k-fold cross-validation for protein subcellular localization problem.

### cmd
```R
Rscript hw3_studentID.R --fold k --input Archaeal_tfpssm.csv --output performance.csv
```
* Perform *k*-fold cross-validation
* % of training, % of calibration, % of testing= *k*-2, 1, 1

The following shows the example of the 5-fold cross validation.
![cross-validation](/images/img2.png)

## Input: Archaeal_tfpssm.csv

This CSV doesn't contain a header. The information of columns as below:

V2: labels of proteins

* CP: Cytoplasmic
* CW: Cell Wall
* EC: Extracellular
* IM: Inner membrane

V3 ~ V5602: the gapped-dipeptide features of each protein

### Code for reference

```R
library('rpart')
# read input data
d <- read.csv(<Path to Archaeal_tfpssm.csv>, header = F)
# label to be predicted
levels(d[,2])
head(d[,5600:5603])
# select subset of the data
tmp <- d[c(seq(1,700,25), seq(700,800,5)),]
# model using decision tree
model <- rpart(V2 ~ V3 + V4 + V5600 + V5601 + V5602,
               data=tmp, control=rpart.control(maxdepth=4),
               method="class")
# make confusion matrix tabel
resultframe <- data.frame(truth=tmp$V2,
                          pred=predict(model, type="class"))
(rtab <- table(resultframe)) 
```

## Model

* Any model you want
* Predict V2 value for each protein

## Output: performance.csv

* accuracy = *P*/*N*, average of *k*-fold cross-validation

set|training|validation|test
---|---|---|---
fold1|0.93|0.91|0.88
fold2|0.92|0.91|0.89
fold3|0.94|0.92|0.90
fold4|0.91|0.89|0.87
fold5|0.90|0.92|0.87
ave.|0.92|0.91|0.88

## Score

* 6 testing cmds from 5-fold to 10-fold
```R
Rscript hw3_studentID.R --fold 5 --input Archaeal_tfpssm.csv --output hw4/your_ID/output1.csv
...
Rscript hw3_studentID.R --fold 10 --input Archaeal_tfpssm.csv --output hw4/your_ID/output6.csv
```
Each testing cmd gets 15 points.
**Please do not set input/output in your local path or URL.** 
Otherwise, your code will fail due to fixed path problem.

## Bonus
* Round number to two decimal places: 3 points
* Without “”: 2 points
* Performance Bonus: average testing performance
  * 0.99 ~: 5 points
  * 0.95 ~ 0.99: 4 points
  * 0.90 ~ 0.95: 3 points
  * 0.80 ~ 0.90: 2 points

## References
* Chang, J.-M. M. et al. (2013) [Efficient and interpretable prediction of protein functional classes by correspondence analysis and compact set relations](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075542). *PLoS ONE* 8, e75542.
* Chang J-M, Su EC-Y, Lo A, Chiu H-S, Sung T-Y, & Hsu W-L (2008) [PSLDoc: Protein subcellular localization prediction based on gapped-dipeptides and probabilistic latent semantic analysis](https://onlinelibrary.wiley.com/doi/full/10.1002/prot.21944). *Proteins: Structure, Function, and Bioinformatics* 72(2):693-710.
