# RECOV
Meta analysis using explicit covariance matrix>


This meta-analysis depends on the traditional random effects (RE) meta-analysis. But, instead of using the tradional identity matrix, it uses a covariance matrix to model the fact that the same SNP can affect the same gene in similar tissues. 

Download the code. 
```
fitLlh1Gene4github.R 
```

Download an example input. 
```
ENSG00000196549.6.beta.txt 
```
There are 44 tissues in the input. The first column is the SNP name. The second column is the SNP effect in tissue 1, and the third column is the variance of this SNP effect in tissue 1. The fourth column is the SNP effect in tissue 2, and the fifth column is the variance of this SNP effect in tissue 2, and so forth. 

Here is an example to run the code 

```
R CMD BATCH --no-save --no-restore "--args /u/home/d/datduong/ENSG00000196549.6.beta.txt /u/home/d/datduong/ENSG00000196549.6.output.txt 10" /u/home/d/datduong/fitLlh1Gene4github.R /u/home/d/datduong/log.txt
```

The entry "10" in the args above splits the SNPs into 10 segments. When fitting the likelihood at a SNP in segment 1, we use the SNPs in the other 9 segments to estimate a rough form of the covariance in the likelihood. 

