# RECOV
Meta analysis using explicit covariance matrix


This meta-analysis depends on the existing random effects (RE) meta-analysis. Instead of using the tradional identity matrix, it uses a covariance matrix to model the fact that the same SNP can affect the same gene in similar tissues. 

Here is an example to run the code 

```
R CMD BATCH --no-save --no-restore "--args /u/home/d/datduong/ENSG00000196549.6.beta.txt /u/home/d/datduong/ENSG00000196549.6.output.txt 10" /u/home/d/datduong/fitLlh1Gene4github.R /u/home/d/datduong/log.txt
```
