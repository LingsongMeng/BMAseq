# BMAseq

BMAseq is an R package for dealing with differential gene expression analysis in heteregenous RNA-seq data.

## Introduction

BMAseq package was developed based on a Bayesian Model Averaging model, which offers many variants on differential gene expression analysis using RNA-seq data. Here, we present multiple analytical scenariors for BMAseq, including univariate analysis, multivariate analysis without interaction, multivariate analysis with interactions, and additional analysis with specific model space settings. More details can be found un the [vignette file](https://github.com/LingsongMeng/BMAseq/blob/master/vignettes/BMAseq_Vignette_V9.Rmd).

## Installation

```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("LingsongMeng/BMAseq")
```

## Basic Usage
### Load sample data
```r
library(BMAseq)

# read sample data
dat.pheno <- dget("/Users/fred/Table/BayesResearch/sample_dat_pheno")
dat.expr <- dget("/Users/fred/Table/BayesResearch/sample_dat_expr")
```

### Univariate analysis 
```r
# perform univariate analysis
vars <- c("BMI", "SEX", "MHHTN", "MHT2D", "MHCVD", "MHHRTATT") 

output.uni <- BMAseq.uni(dat.expr.counts = dat.expr, 
                         dat.pheno = dat.pheno, 
                         var.pool = vars, 
                         cut.BF = 1, 
                         cut.FDR = 0.5)

# return the number of the differentially expressed genes associated with each variable of interest
output.uni$nDEG

# return the differentially expressed genes associated with each variable of interest
output.uni$DEG$BMI
output.uni$DEG$SEX
output.uni$DEG$MHHTN
output.uni$DEG$MHT2D
output.uni$DEG$MHCVD
output.uni$DEG$MHHRTATT
```


### Multivariate analysis without interaction 
```r
# perform multivariate analysis without interaction
var.pool <- c("BMI", "SEX", "MHHTN", "MHT2D") 

output.multi <- BMAseq.multi(dat.expr.counts = dat.expr, 
                             dat.pheno = dat.pheno, 
                             var.pool = var.pool, 
                             max.nvar = 4, 
                             interaction = NULL, 
                             cut.BF = 1, 
                             cut.FDR = 0.05)

# return a summary table of the number of the identified differentially expressed gene associated with the main effect of each variable
output.multi$summary.nDEG

# return a summary table of the number of the identified differentially expressed genes associated with the joint main effects of variables
output.multi$summary.nDEG.JointMain

# return the differentially expressed genes associated with the main effect of each variable, and the best model used to identify each differentially expressed gene  
output.multi$DEG.bestmodel   
```


### Multivariate Analysis with Interaction
```r
# perform multivariate analysis with interaction
var.pool <- c("BMI", "SEX", "MHHTN", "MHT2D") 

interaction <- "BMI&SEX"  

output.multi.int1 <- BMAseq.multi.postprob(dat.expr.counts = dat.expr, 
                                           dat.pheno = dat.pheno, 
                                           var.pool = var.pool, 
                                           max.nvar = 4, 
                                           interaction = interaction, 
                                           cut.BF = 1)

output.multi.int2 <- BMAseq.multi.DEG(postprob.output = output.multi.int1, 
                                      cut.FDR = 0.05)

# return a summary table of the number of the identified differentially expressed gene associated with each variable
output.multi.int2$summary.nDEG

# return the differentially expressed genes associated with the main, interaction, main or interaction effect of each variable, and the best model used to identify each differentially expressed gene
output.multi.int2$DEG.bestmodel
```

## Advanced Usage
### Additional Analysis with Specific Model Space Settings
```r
# build a model sapce and calculate posterior model probability
MS_default <- Modelspace(dat.pheno = dat.pheno, 
                         var.pool = var.pool, 
                         max.nvar = 4, 
                         interaction = interaction)

output.multi.int1 <- BMAseq.multi.postprob.MSout(dat.expr.counts = dat.expr, 
                                                 dat.pheno = MS_default$dat.pheno.new,
                                                 model.space = MS_default$model.space, 
                                                 var.pool = MS_default$var.pool,
                                                 interaction = MS_default$interaction,
                                                 cut.BF = 1)

# identify the differentially expressed genes associated with the additional specific inclusive models
# 1. BMI (all the models containing BMI but no BMI related interaction terms)
id1 <- c(2, 4, 18, 20, 34, 36, 50, 52)
# 2. SEX (all the models containing SEX but no SEX related interaction terms)
id2 <- c(3, 4, 19, 20, 35, 36, 51, 52)
# 3. BMIlow.SEXfemale (all the models with this interaction term but not BMI or SEX as main effects)
id3 <- c(5, 21, 37, 53)
# 4. BMIlow.SEXmale (all the models with this interaction term but not BMI or SEX as main effects)
id4 <- c(6, 22, 38, 54)
# 5. BMIhigh.SEXfemale (all the models with this interaction term but not BMI or SEX as main effects)
id5 <- c(7, 14, 23, 30, 39, 46, 55, 62)
# 6. BMIhigh.SEXmale (all the models with this interaction term but not BMI or SEX as main effects)
id6 <- c(8, 15, 24, 31, 40, 47, 56, 63)
# 7. BMIlow.SEXfemale.BMIhigh.SEXmale (all the models with this interaction term but not BMI or SEX as main effects)
id7 <- c(9, 14, 15, 25, 30, 31, 41, 46, 47, 57, 62, 63)
# 8. BMI + BMIlow.SEXmale
id8 <- c(10, 26, 42, 58)
# 9. BMI + BMIhigh.SEXmale
id9 <- c(12, 28, 44, 60)
# 10. SEX + BMIhigh.SEXfemale
id10 <- c(11, 27, 43, 59)
# 11. SEX + BMIhigh.SEXmale
id11 <- c(13, 29, 45, 61)
# 12. BMIlow.SEXfemale.BMIhigh.SEXmale + BMIhigh.SEXfemale
id12 <- c(14, 30, 46, 62)
# 13. BMIlow.SEXfemale.BMIhigh.SEXmale + BMIhigh.SEXmale
id13 <- c(15, 31, 47, 63)
# 14. BMI + SEX
id14 <- c(4, 20, 36, 52)
# 15. BMI + SEX + BMIhigh.SEXmale
id15 <- c(16, 32, 48, 64)

id.add.list <- list(id1, id2, id3, id4, id5, id6, id7, id8, id9, id10, id11, id12, id13, id14, id15)
name <- c("BMI.Main.noInt", "SEX.Main.noInt", "BMIlowSEXfemale.noMain", "BMIlowSEXmale.noMain",  "BMIhighSEXfemale.noMain", "BMIhighSEXmale.noMain", "BMIlowSEXfemaleBMIhighSEXmale.noMain", "BMI.plus.BMIlowSEXmale", "BMI.plus.BMIhighSEXmale", "SEX.plus.BMIhighSEXfemale", "SEX.plus.BMIhighSEXmale", "BMIlowSEXfemaleBMIhighSEXmale.plus.BMIhighSEXfemale", "BMIlowSEXfemaleBMIhighSEXmale.plus.BMIhighSEXmale", "BMI.plus.SEX", "BMI.plus.SEX.plus.BMIhighSEXmale")
names(id.add.list) <- name


output.multi.int2.add <- BMAseq.multi.DEG(postprob.output = output.multi.int1,
                                          ind.incl.add = id.add.list, 
                                          cut.FDR = 0.05)

# return a summary table of the number of the identified differentially expressed genes associated with the additionally specified models
summary.nDEG.add
```

# Getting help

If you have any inquiries regarding our BMAseq, don't hesitate to reach out to us at lingsongmeng2019@gmail.com and anni.liu.001@gmail.com.

# Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.
