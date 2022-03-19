# Bootstrapping the t-test
Bootstrap_di <- function (Losses,boot_samples,testtype, N){
  # Constructs that has the same name in both models
  same_constructs_both_models <- intersect(names(Losses$LossM1_sepLV),names(Losses$LossM2_sepLV))
  # Combine losses in one dataframe
  all_loss <- all_loss_in_df(Losses, same_constructs_both_models)
  # All t-test in one df
  all_Ttest <- all_Ttest_in_df(all_loss, same_constructs_both_models, testtype)
  # All original average difference in losses
  all_Dbar <- all_Dbar_in_df(all_loss, same_constructs_both_models)
  # All model differences in loss functions under the null
  all_D_0 <- all_D_0_in_df(all_loss, all_Dbar, same_constructs_both_models)
  # All model differences in loss functions
  all_D <- all_D_in_df(all_loss, same_constructs_both_models)
  #Allocating memory to bootstrap
  BootSample <- matrix(0,ncol=ncol(all_loss),nrow=nrow(all_D))
  BootDbar <- matrix(0,ncol=(1+length(same_constructs_both_models)),nrow=boot_samples,
                     dimnames = list(NULL,c("Overall",same_constructs_both_models)))
  tStat <- matrix(0,ncol=(1+length(same_constructs_both_models)),nrow=boot_samples,
                  dimnames = list(NULL,c("Overall",same_constructs_both_models)))
  # Bootstrapping
  for (b in 1:boot_samples) {
    # Random draw of indices with replacement
    # TODO: We could make only one draw - better aligning the results. However, this would make a small disrepancy
    # with respect to the "predictive model assessment paper".
    random_indices1 <- sample((1:nrow(all_D)), nrow(all_D), replace=TRUE)
    random_indices2 <- sample((1:nrow(all_D_0)), nrow(all_D_0), replace=TRUE)
    # Draw from original losses
    BootSample <-  all_loss[random_indices1,]
    # t-test on bootstrapped sample, using original D_bar to test under null
    tStat[b,]<-unlist(all_boot_Ttest_in_df(BootSample, all_D, same_constructs_both_models, testtype))
    # Random draw on loss differences under the null hypothesis
    BootDbar[b,] <- colMeans(all_D_0[random_indices2,])
  }
  # Sort bootstrapped t-statistic for overall model and each construct
  SorttStat <- tStat
  for (c in colnames(tStat)) {
    SorttStat[,c] <- sort(tStat[,c], decreasing = FALSE)
  }
  # Sort bootstrapped D-bar
  SortBootDbar <- BootDbar
  for (c in colnames(BootDbar)) {
    SortBootDbar[,c] <- sort(BootDbar[,c], decreasing = FALSE)
  }
  # Bootstrap variance on D-bar for t-test
  all_std <- apply(BootDbar, MARGIN = 2, function(x) sqrt(var(x)) )
  # t-statistic using bootstrapped variance of D-bar
  all_tstat_boot_Var <- all_Dbar/all_std
  colnames(all_tstat_boot_Var) <- gsub("OrgDbar", "tstat_boot_Var", colnames(all_tstat_boot_Var))
  # Calculating p-values
  if (testtype=="two.sided") {
    # Fraction of bootstrapped t-test more extreme than original t-test
    p.value_perc_Ttest <- mapply(FUN = function(x,y) {
      (sum(SorttStat[,x]>abs(all_Ttest[,y]))+sum(SorttStat[,x]<=(-abs(all_Ttest[,y]))))/boot_samples
    },
    colnames(SorttStat), colnames(all_Ttest))
    # Fraction of bootstrapped D-bar more extreme than original D-bar
    p.value_perc_D <- mapply(FUN = function(x,y) {
      (sum(SortBootDbar[,x]>abs(all_Dbar[,y]))+sum(SortBootDbar[,x]<=(-abs(all_Dbar[,y]))))/boot_samples
    },
    colnames(SortBootDbar), colnames(all_Dbar))
    # Evaluating t-stat using bootstrapped variance against t-distribution with N-1 df
    p.value_var_ttest <- unlist(lapply(all_tstat_boot_Var,
                                       function(x) {2*pt(-abs(x),(N-1), lower.tail = TRUE)}))
  }

  if (testtype=="greater") {
    # Fraction of bootstrapped t-test larger than original t-test
    p.value_perc_Ttest <- mapply(FUN = function(x,y) {
      if (!any(which(SorttStat[,x]>all_Ttest[,y]))) {
        0
      }else{
        1-(head(which(SorttStat[,x]>all_Ttest[,y]),1)-1)/(boot_samples+1)
      }
    },
    colnames(SorttStat), colnames(all_Ttest))
    # Fraction of bootstrapped D-bar larger than original D-bar
    p.value_perc_D <- mapply(FUN = function(x,y) {
      if (!any(SortBootDbar[,x]>all_Dbar[,y])) {
        0
      }else{
        1-(head(which(SortBootDbar[,x]>all_Dbar[,y]),1)-1)/(boot_samples+1)
      }
    },
    colnames(SortBootDbar), colnames(all_Dbar))
    # Evaluating t-stat using bootstrapped variance against t-distribution with N-1 df
    p.value_var_ttest <- unlist(lapply(all_tstat_boot_Var,
                                       function(x) {pt(x,(N-1), lower.tail = FALSE)}))

    # If none of the bootstrapped t-tests is larger than the original
    # for (a in names(p.value_perc_Ttest)) {
    #   if (length(which(SorttStat[,a]>all_Ttest[,paste("OrgTtest", a, sep = "_")]))==0){
    #     p.value_perc_Ttest[a]=0
    #     p.value_perc_D[a]=0
    #   }
    #
    #
    # }


  }
  # Naming and classes uniform across results
  all_Ttest <- unlist(all_Ttest)
  names(all_Ttest) <- names(p.value_perc_D)

  all_tstat_boot_Var <- unlist(all_tstat_boot_Var)
  names(all_tstat_boot_Var) <- names(p.value_perc_D)

  names(p.value_var_ttest) <- names(p.value_perc_D)

  # List of all results
  Results <- list("OrgTtest"=all_Ttest, "p.value.perc.t"=p.value_perc_Ttest,
               "t.stat.b.v" = all_tstat_boot_Var, "p.value.var.ttest" = p.value_var_ttest,
               "p.value.perc.D" = p.value_perc_D)
  return(Results)
}

## Helpers for bootstrapping_di
all_loss_in_df <- function(Losses, same_constructs_both_models){
  # Initiate df with overall losses
  all_loss <- data.frame(LossM1 = Losses$LossM1, LossM2 = Losses$LossM2)
  # Insert losses from each model on same constructs
  for (c in same_constructs_both_models) {
    all_loss[,paste("LossM1", c, sep = "_")] <- Losses$LossM1_sepLV[[c]]
    all_loss[,paste("LossM2", c, sep = "_")] <- Losses$LossM2_sepLV[[c]]
  }
  return(all_loss)
}

all_Ttest_in_df <- function(all_loss, same_constructs_both_models, testtype){
  # T-test on overall model
  OrgTtest_Overall<-t.test(all_loss$LossM2,all_loss$LossM1,alternative = testtype, paired=TRUE)$statistic
  # Initiate df with overall t-test
  all_Ttest <- data.frame(OrgTtest_Overall = OrgTtest_Overall)
  # t-test on each construct (using loss on each construct from seperate models)
  for (c in same_constructs_both_models) {
    all_Ttest[,paste("OrgTtest", c, sep = "_")] <- t.test(all_loss[,paste("LossM2", c, sep = "_")],
                                                          all_loss[,paste("LossM1", c, sep = "_")],
                                                          alternative = testtype, paired=TRUE)$statistic
  }
  return(all_Ttest)
}

all_Dbar_in_df <- function(all_loss, same_constructs_both_models){
  # Dbar on overall model
  OrgDbar_overall <- mean(all_loss$LossM2-all_loss$LossM1)
  # Initiate df with overall Dbar
  all_OrgDbar <- data.frame(OrgDbar_overall = OrgDbar_overall)
  # Dbar on each construct (using loss on each construct from seperate models)
  for (c in same_constructs_both_models) {
    all_OrgDbar[,paste("OrgDbar", c, sep = "_")] <- mean(all_loss[,paste("LossM2", c, sep = "_")] -
                                                           all_loss[,paste("LossM1", c, sep = "_")])
  }
  return(all_OrgDbar)
}

all_D_0_in_df <- function(all_loss, all_Dbar, same_constructs_both_models){
  # Dbar on overall model
  D_0_overall <- all_loss$LossM2 - all_loss$LossM1 - all_Dbar$OrgDbar_overall
  # Initiate df with overall Dbar
  all_D_0 <- data.frame(D_0_overall = D_0_overall)
  # Dbar on each construct (using loss on each construct from seperate models)
  for (c in same_constructs_both_models) {
    all_D_0[,paste("D_0", c, sep = "_")] <- all_loss[,paste("LossM2", c, sep = "_")] -
      all_loss[,paste("LossM1", c, sep = "_")] -
      all_Dbar[,paste("OrgDbar", c, sep = "_")]
  }
  return(all_D_0)
}

all_D_in_df <- function(all_loss, same_constructs_both_models){
  # Dbar on overall model
  D_overall <- all_loss$LossM2 - all_loss$LossM1
  # Initiate df with overall D_i
  all_D <- data.frame(Overall = D_overall)
  # D on each construct (using loss on each construct from seperate models)
  for (c in same_constructs_both_models) {
    all_D[,c] <- all_loss[,paste("LossM2", c, sep = "_")] -
      all_loss[,paste("LossM1", c, sep = "_")]
  }
  return(all_D)
}

all_boot_Ttest_in_df <- function(all_loss, all_D, same_constructs_both_models, testtype){
  # T-test on overall model
  OrgTtest_Overall<-t.test(all_loss$LossM2,all_loss$LossM1, mu = mean(all_D[,"Overall"]),
                           alternative = testtype, paired=TRUE)$statistic
  # Initiate df with overall t-test
  all_Ttest <- data.frame(OrgTtest_Overall = OrgTtest_Overall)
  # t-test on each construct (using loss on each construct from seperate models)
  for (c in same_constructs_both_models) {
    all_Ttest[,paste("OrgTtest", c, sep = "_")] <- t.test(all_loss[,paste("LossM2", c, sep = "_")],
                                                          all_loss[,paste("LossM1", c, sep = "_")],
                                                          mu = mean(all_D[,c]),
                                                          alternative = testtype, paired=TRUE)$statistic
  }
  return(all_Ttest)
}

# Bootstrap on MV
# Bootstrap_MV <- function(mv.org.ttest, Di, D.bar, N, testtype, MV,CVFolds,InnerSpecM1,ReflecSpecM1,FormSpecM1,InnerSpecM2,ReflecSpecM2,FormSpecM2, boot_samples){
#   b.t.stat <- rep(0,boot_samples)
#   b.D.bar <- rep(0,boot_samples)
#   # Run bootstrapping
#   for (b in 1:boot_samples) {
#     b.samp<-sample(1:N,N,replace = TRUE)
#     b.MV<-MV[b.samp,]
#     # Use trycatch, to count proportion of times with non-convergence
#     b.Losses <- tryCatch(LossFun(N,b.MV,CVFolds,InnerSpecM1,ReflecSpecM1,FormSpecM1,InnerSpecM2,ReflecSpecM2,FormSpecM2),
#                          warning = function(w) {NA},
#                          error = function(e) {NA})
#     # Save as NA if non-convergence, else save bootstrapped calculated values
#     if (is.na(b.Losses[1])) {
#       b.D.bar[b] <- NA
#       b.t.stat[b]<- NA
#     } else {
#       b.D.bar[b] <- mean(b.Losses$LossM2 - b.Losses$LossM1) - D.bar
#       b.t.stat[b]<-t.test(b.Losses$LossM2,b.Losses$LossM1,mu=D.bar,alternative=testtype, paired=TRUE)$statistic
#     }
#   }
#   # Proportion of non-convergent bootstrap runs
#   prop.non.conv <- sum(is.na(b.D.bar))/length(b.D.bar)
#   # Number of bootstrap runs that did converge
#   boot_samples<-length(b.D.bar[!is.na(b.D.bar)])
#   # Discarding non-convergent bootstrap runs
#   b.D.bar<-b.D.bar[!is.na(b.D.bar)]
#   b.t.stat<-b.t.stat[!is.na(b.t.stat)]
#
#   # Sorted bootstrapped t-statistics
#   t.sort <- sort(b.t.stat, decreasing = FALSE)
#   # Bootstrapped variance on D_bar
#   t.stat.b.v <- D.bar/sqrt(var(b.D.bar))
#   # Sorted bootstrapped D_bar
#   D.bar.sort <- sort(b.D.bar, decreasing = FALSE)
#   # P-values for two-sided test
#   if (testtype=="two.sided") {
#     # Percentile-t method, p-value
#     p.value.perc.t<-(sum(t.sort>abs(mv.org.ttest))+sum(t.sort<=(-abs(mv.org.ttest))))/boot_samples
#     # Botstrapped variance of D.bar, p-value
#     p.value.var.ttest<-2*pt(-abs(t.stat.b.v),(N-1), lower.tail = TRUE)
#     # percentile D.bar, p-value
#     p.value.perc.D<-(sum(D.bar.sort>abs(D.bar))+sum(D.bar.sort<=(-abs(D.bar))))/boot_samples
#   }
#   if (testtype=="greater") {
#     # Percentile-t method, p-value
#     p.value.perc.t <- 1-(head(which(t.sort>mv.org.ttest),1)-1)/(boot_samples+1)
#     if (length(which(t.sort>mv.org.ttest))==0){
#       p.value.perc.t=0
#     }
#     # Botstrapped variance of D.bar, p-value
#     p.value.var.ttest <- pt(t.stat.b.v,(N-1), lower.tail = FALSE)
#     # percentile D.bar, p-value
#     p.value.perc.D <- 1-(head(which(D.bar.sort>D.bar),1)-1)/(boot_samples+1)
#     if (length(which(D.bar.sort>D.bar))==0){
#       p.value.perc.D=0
#     }
#   }
#   return(c("p.value.perc.t" = p.value.perc.t, "p.value.var.ttest" = p.value.var.ttest, "p.value.perc.D" = p.value.perc.D,
#            "t.stat.b.v" = t.stat.b.v, "prop.non.conv" = prop.non.conv))
# }
