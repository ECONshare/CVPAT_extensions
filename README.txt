###########################################################
## CITATION: Sharma, P.N., Liengaard, B., Hair, J., Sarstedt, M., and Ringle, C. "Predictive model assessment and selection in composite-based modeling using PLS-SEM: extensions and guidelines for using CVPAT" (European Journal of Marketing, 2022).
###########################################################

###########################################################
## Input to CVPAT function
###########################################################
 
mv: 
Data of the manifest variables (indicators) used for the partial least squares (PLS) path model estimation. The variable names must be the same as specified in model1 and model2 (see below).

cv_folds:
Sets the number of cross-validation folds.

model1: 
either:
1) The PLS model specifications as required by the cSEM package (see also the CVPAT_example.R file) or
2) "IA" for indicator average benchmark or
3) "LM" for the linear model benchmark

Model2: 
either:
1) The PLS model specifications as required by the cSEM package (see also the CVPAT_example.R file) or
2) "IA" for indicator average benchmark or
3) "LM" for the linear model benchmark

hypothesis: 
Specifies the hypothesis to be tested. If hypothesis = "M1_better_out_of_sample_than_M2", then CVPAT will test whether M1 has a significantly better performance out-of-sample performance than M2. 
If Hypothesis = "M1!=M2", then CVPAT will test if the out-of-sample performance between M1 and M2 differs significantly from each other.

compare_on_constructs:
A character vector with the names of the reflectively measured endogenous latent variables for which the CVPAT is performed.

boot_samples:
Sets the number of bootstrap samples.

seed: 
Sets the seed of the random number generator.



###########################################################
## Output from cvpat function
###########################################################
The output from the CVPAT function is a list. The list has following elements:
- res
- conf.int
- losses

----------
The content of each list element is:
----------

res:
A matrix with column names for endogenous latent variables. If a set of constructs is given in the "compare_on_constructs" argument then the results
for these constructs combined will be given in the last column (the construct names will be seperated by "_").
- avg_loss_M1: The average loss for model 1
- avg_loss_M2: The average loss for model 2
- Diff (M2 - M1): The difference betweeen the average loss in model 1 and model 2.
- non.boot.t.stat: Non bootstrapped t-statistic
- non.boot.p.val: Non bootstrapped p-value
- boot.var.t.stat: t-statistics using bootstrapped variance
- p.value.perc.t: p-value calculated from the bootstrapped percentiles of the t-statistics.
- p.value.b.v.t: p-value calculated from the original t-statistic but replacing the variance of D_bar with the bootstrapped variance of D_bar
- p.value.perc.D: p-value calculated from the bootstrapped percentiles of D_bar.

conf.int: 
The non-bootstrapped confidence interval

losses: 
List of losses for various subparts of the model

