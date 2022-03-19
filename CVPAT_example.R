#IMPORTANT: SET YOUR PATH TO THE WORKING DIRECTORY
#WHERE THE files "corp_rep.rda", "cvpat.R", "loss_function.R" and "bootstrap_functions.R" ARE LOCATED
setwd("C:/CVPAT")


# Install required packages --------------------------------------------------
#install.packages("cSEM")
#install.packages("cvTools")


# Load required packages --------------------------------------------------
library(cSEM)
library(cvTools)

# Load data ---------------------------------------------------------------
load(file="corp_rep.rda")
# Source scripts ----------------------------------------------------------
source("cvpat.R")
source("loss_function.R")
source("bootstrap_functions.R")

# Specify all models ---------------------------------------------
## Established Model (EM)
EM <- "
# Structural model
COMP ~ QUAL + PERF + CSOR + ATTR
LIKE ~ QUAL + PERF + CSOR + ATTR
CUSA ~ COMP + LIKE
CUSL ~ COMP + CUSA + LIKE

# Composite model
QUAL <~ qual_1 + qual_2 + qual_3 + qual_4 + qual_5 + qual_6 + qual_7 + qual_8
PERF <~ perf_1 + perf_2 + perf_3 + perf_4 + perf_5
CSOR <~ csor_1 + csor_2 + csor_3 + csor_4 + csor_5
ATTR <~ attr_1 + attr_2 + attr_3

# Reflective measurement model
COMP  =~ comp_1 + comp_2 + comp_3
LIKE  =~ like_1 + like_2 + like_3
CUSA  =~ cusa
CUSL  =~ cusl_1 + cusl_2 + cusl_3
"

## Alternative Model 1 (AM1)
AM1 <- "
# Structural model
COMP ~ QUAL + PERF
LIKE ~ CSOR + ATTR
CUSA ~ COMP + LIKE
CUSL ~ COMP + CUSA + LIKE

# Composite model
QUAL <~ qual_1 + qual_2 + qual_3 + qual_4 + qual_5 + qual_6 + qual_7 + qual_8
PERF <~ perf_1 + perf_2 + perf_3 + perf_4 + perf_5
CSOR <~ csor_1 + csor_2 + csor_3 + csor_4 + csor_5
ATTR <~ attr_1 + attr_2 + attr_3

# Reflective measurement model
COMP  =~ comp_1 + comp_2 + comp_3
LIKE  =~ like_1 + like_2 + like_3
CUSA  =~ cusa
CUSL  =~ cusl_1 + cusl_2 + cusl_3
"
## Alternative Model 2 (AM2)
AM2 <- "
# Structural model
COMP ~ QUAL + PERF + CSOR + ATTR
LIKE ~ QUAL + PERF + CSOR + ATTR
CUSA ~ COMP + LIKE + QUAL + PERF + CSOR + ATTR
CUSL ~ COMP + CUSA + LIKE + QUAL + PERF + CSOR + ATTR


# Composite model
QUAL <~ qual_1 + qual_2 + qual_3 + qual_4 + qual_5 + qual_6 + qual_7 + qual_8
PERF <~ perf_1 + perf_2 + perf_3 + perf_4 + perf_5
CSOR <~ csor_1 + csor_2 + csor_3 + csor_4 + csor_5
ATTR <~ attr_1 + attr_2 + attr_3

# Reflective measurement model
COMP  =~ comp_1 + comp_2 + comp_3
LIKE  =~ like_1 + like_2 + like_3
CUSA  =~ cusa
CUSL  =~ cusl_1 + cusl_2 + cusl_3
"

# CVPAT Benchmark -----------------------------------------------------
## Focus on CUSA and CUSL constructs

# EM and Indicator Average benchmark
res_CVPAT_IA <- cvpat(mv = corp_rep, cv_folds = 10,
                   model1 = EM, model2 = "IA",
                   hypothesis = "M1_better_out_of_sample_than_M2",
                   compare_on_constructs = c("CUSA", "CUSL"),
                   boot_samples= 2000,
                   rescale_pred_errors = FALSE,
                   seed=42)
# Inspect CVPAT results
res_CVPAT_IA$res[c("avg_loss_M1", "avg_loss_M2", "Diff (M2 - M1)", "p.value.perc.t"),
                 "CUSA_CUSL"]

# EM and Indicator Average benchmark
res_CVPAT_LM <- cvpat(mv = corp_rep, cv_folds = 10,
                   model1 = EM, model2 = "LM",
                   hypothesis = "M1_better_out_of_sample_than_M2",
                   compare_on_constructs = c("CUSA", "CUSL"),
                   boot_samples= 2000,
                   rescale_pred_errors = FALSE,
                   seed=42)

# Inspect CVPAT results
res_CVPAT_LM$res[c("avg_loss_M1", "avg_loss_M2", "Diff (M2 - M1)", "p.value.perc.t"),
                 "CUSA_CUSL"]

# CVPAT compare EM with AM1 -----------------------------------------------------
## Focus on CUSA and CUSL constructs

# EM and Indicator Average benchmark
res_CVPAT_AM1 <- cvpat(mv = corp_rep, cv_folds = 10,
                      model1 = AM1, model2 = EM,
                      hypothesis = "M1_better_out_of_sample_than_M2",
                      compare_on_constructs = c("CUSA", "CUSL"),
                      boot_samples= 2000,
                      rescale_pred_errors = FALSE,
                      seed=42)
# Inspect CVPAT results
res_CVPAT_AM1$res[c("avg_loss_M1", "avg_loss_M2", "Diff (M2 - M1)", "p.value.perc.t"),
                 "CUSA_CUSL"]
# CVPAT compare EM with AM2 -----------------------------------------------------
# EM and Indicator Average benchmark
res_CVPAT_AM2 <- cvpat(mv = corp_rep, cv_folds = 10,
                      model1 = AM2, model2 = EM,
                      hypothesis = "M1_better_out_of_sample_than_M2",
                      compare_on_constructs = c("CUSA", "CUSL"),
                      boot_samples= 2000,
                      rescale_pred_errors = FALSE,
                      seed=42)

# Inspect CVPAT results
res_CVPAT_AM2$res[c("avg_loss_M1", "avg_loss_M2", "Diff (M2 - M1)", "p.value.perc.t"),
                 "CUSA_CUSL"]


