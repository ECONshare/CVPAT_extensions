LossFun <- function(N, mv, cv_folds, model1, model2, compare_on_constructs, rescale_pred_errors) {
  mod_comparison_type <- check_and_classify_fun_input(N = N,
                                                      mv = mv,
                                                      cv_folds = cv_folds,
                                                      model1 = model1,
                                                      model2 = model2)
  # Estimate PLS models
  switch(mod_comparison_type,
         m1_pls_m2_bench = {
           # Estimate pls
           model1_pls <- est_pls_fun(mv = mv,
                                     model = model1)
           # Memory for prediction errors
           PredErrorsM1_MV <- PredErrors_IA <- PredErrors_LM <- mv[,model1_pls$ReflEndogMV]*0
         },
         m1_bench_m2_pls = {
           # Estimate pls
           model2_pls <- est_pls_fun(mv = mv,
                                     model = model2)
           # Memory for prediction errors
           PredErrorsM2_MV <- PredErrors_IA <- PredErrors_LM <- mv[,model2_pls$ReflEndogMV]*0
         },
         m1_pls_m2_pls = {
           # Estimate pls
           model1_pls <- est_pls_fun(mv = mv,
                                     model = model1)
           model2_pls <- est_pls_fun(mv = mv,
                                     model = model2)
           # Check if same constructs is compared between models
           common_endog_LV <- intersect(model1_pls$ReflEndogLV, model2_pls$ReflEndogLV)
           actual_compare <- intersect(common_endog_LV, compare_on_constructs)
           if(!all(compare_on_constructs %in% actual_compare)){
             # Only compare constructs that appear in both models
             compare_on_constructs <- intersect(common_endog_LV, compare_on_constructs)
             message("
             The constructs you specified in the compare_on_construct argument is either \n
             1) Not endogenous or reflective in both models \n
             2) Misspelled
             Please check your input to the compare_on_construct argument11.")

             if (length(compare_on_constructs)==0) {
               compare_on_constructs <- NULL
               warning(paste("No construct can be compared between the models"))
             }
             warning(paste("Only construct",paste(compare_on_constructs, collapse = " "),
                           "can be compared between the models"))

           }


           # Memory for prediction errors
           PredErrorsM1_MV <- matrix(0, nrow = nrow(mv), ncol = length(model1_pls$ReflEndogMV),
                                     dimnames = list(NULL, model1_pls$ReflEndogMV))
           PredErrorsM2_MV <- matrix(0, nrow = nrow(mv), ncol = length(model2_pls$ReflEndogMV),
                                     dimnames = list(NULL, model2_pls$ReflEndogMV))
         }
  )

  ## Split data into folds

  folds<-cvTools::cvFolds(N,cv_folds,1)
  folds<-cbind(folds$which,folds$subsets)
  for (i in 1:cv_folds) {
    holdoutIndexes <- folds[folds[,1]==i,2]
    FoldSize <- length(holdoutIndexes)
    TrainMV     <- mv[-holdoutIndexes,]
    ## Standardize training data
    Train_mean <- colMeans(TrainMV)
    Train_std <- sapply(colnames(TrainMV), FUN = function(x){sqrt(var(TrainMV[,x]))})
    TrainMV <- sapply(colnames(TrainMV), FUN = function(x){(TrainMV[,x] - Train_mean[x])/Train_std[x]})
    # Standardized holdout data
    HoldoutMV <- sapply(colnames(mv[holdoutIndexes,]),
                        FUN = function(x){(mv[holdoutIndexes,x] - Train_mean[x])/Train_std[x]},
                        USE.NAMES = FALSE)
    # Naming and ensure that we have a matrix class in the leave-one-out case
    HoldoutMV <- matrix(HoldoutMV,
                        ncol = ncol(mv),
                        dimnames = list(NULL, colnames(mv)))
    # Prediction errors
    switch (mod_comparison_type,
            m1_bench_m2_pls = {
              out_samp_pls_m2 <- out_of_sample_pls(model = model2,
                                                   HoldoutMV =  HoldoutMV,
                                                   FoldSize = FoldSize,
                                                   rescale_pred_errors = rescale_pred_errors,
                                                   TrainMV = TrainMV,
                                                   Train_std = Train_std,
                                                   Train_mean = Train_mean)
              PredErrorsM2_MV[holdoutIndexes,] <- out_samp_pls_m2$PredErrors
              out_samp_bench <- out_of_sample_bench(out_samp_pls = out_samp_pls_m2,
                                                    HoldoutMV =  HoldoutMV,
                                                    FoldSize = FoldSize,
                                                    rescale_pred_errors = rescale_pred_errors,
                                                    TrainMV = TrainMV,
                                                    Train_std = Train_std,
                                                    Train_mean = Train_mean)
              PredErrors_IA[holdoutIndexes,] <- out_samp_bench$PredErrors_IA
              PredErrors_LM[holdoutIndexes,] <- out_samp_bench$PredErrors_LM

              if (model1 == "IA") {
                PredErrorsM1_MV <- PredErrors_IA
              }
              if (model1 == "LM") {
                PredErrorsM1_MV <- PredErrors_LM
              }
            },
            m1_pls_m2_bench = {
              out_samp_pls_m1 <- out_of_sample_pls(model = model1,
                                                   HoldoutMV =  HoldoutMV,
                                                   FoldSize = FoldSize,
                                                   rescale_pred_errors = rescale_pred_errors,
                                                   TrainMV = TrainMV,
                                                   Train_std = Train_std,
                                                   Train_mean = Train_mean)
              PredErrorsM1_MV[holdoutIndexes,] <- out_samp_pls_m1$PredErrors
              out_samp_bench <- out_of_sample_bench(out_samp_pls = out_samp_pls_m1,
                                                    HoldoutMV =  HoldoutMV,
                                                    FoldSize = FoldSize,
                                                    rescale_pred_errors = rescale_pred_errors,
                                                    TrainMV = TrainMV,
                                                    Train_std = Train_std,
                                                    Train_mean = Train_mean)
              PredErrors_IA[holdoutIndexes,] <- out_samp_bench$PredErrors_IA
              PredErrors_LM[holdoutIndexes,] <- out_samp_bench$PredErrors_LM

              if (model2 == "IA") {
                PredErrorsM2_MV <- PredErrors_IA
              }
              if (model2 == "LM") {
                PredErrorsM2_MV <- PredErrors_LM
              }
              },
            m1_pls_m2_pls = {
              out_samp_pls_m1 <- out_of_sample_pls(model = model1,
                                                   HoldoutMV =  HoldoutMV,
                                                   FoldSize = FoldSize,
                                                   rescale_pred_errors = rescale_pred_errors,
                                                   TrainMV = TrainMV,
                                                   Train_std = Train_std,
                                                   Train_mean = Train_mean)
              PredErrorsM1_MV[holdoutIndexes,] <- out_samp_pls_m1$PredErrors
              out_samp_pls_m2 <- out_of_sample_pls(model = model2,
                                                   HoldoutMV =  HoldoutMV,
                                                   FoldSize = FoldSize,
                                                   rescale_pred_errors = rescale_pred_errors,
                                                   TrainMV = TrainMV,
                                                   Train_std = Train_std,
                                                   Train_mean = Train_mean)
              PredErrorsM2_MV[holdoutIndexes,] <- out_samp_pls_m2$PredErrors
              }
    )
  }

  # Losses
  switch (mod_comparison_type,
          m1_bench_m2_pls = {
            LossM2 <- rowMeans(PredErrorsM2_MV^2)
            LossM2_sepLV <- loss_each_lv(model = model2_pls,
                                         PredErrors = PredErrorsM2_MV,
                                         compare_on_constructs = compare_on_constructs)
            Loss_IA <- rowMeans(PredErrors_IA^2)
            Loss_IA_sepLV <- loss_each_lv(model = model2_pls,
                                          PredErrors = PredErrors_IA,
                                          compare_on_constructs = compare_on_constructs)
            Loss_LM <- rowMeans(PredErrors_LM^2)
            Loss_LM_sepLV <- loss_each_lv(model = model2_pls,
                                          PredErrors = PredErrors_LM,
                                          compare_on_constructs = compare_on_constructs)
            if (model1 == "IA") {
              LossM1 <- Loss_IA
              LossM1_sepLV <- Loss_IA_sepLV
            }
            if (model1 == "LM") {
              LossM1 <- Loss_LM
              LossM1_sepLV <- Loss_LM_sepLV
            }
          },
          m1_pls_m2_bench = {
            LossM1 <- rowMeans(PredErrorsM1_MV^2)
            LossM1_sepLV <- loss_each_lv(model = model1_pls,
                                         PredErrors = PredErrorsM1_MV,
                                         compare_on_constructs = compare_on_constructs)
            Loss_IA <- rowMeans(PredErrors_IA^2)
            Loss_IA_sepLV <- loss_each_lv(model = model1_pls,
                                          PredErrors = PredErrors_IA,
                                          compare_on_constructs = compare_on_constructs)
            Loss_LM <- rowMeans(PredErrors_LM^2)
            Loss_LM_sepLV <- loss_each_lv(model = model1_pls,
                                          PredErrors = PredErrors_LM,
                                          compare_on_constructs = compare_on_constructs)
            if (model2 == "IA") {
              LossM2 <- Loss_IA
              LossM2_sepLV <- Loss_IA_sepLV
            }
            if (model2 == "LM") {
              LossM2 <- Loss_LM
              LossM2_sepLV <- Loss_LM_sepLV
            }
          },
          m1_pls_m2_pls = {
            LossM1 <- rowMeans(PredErrorsM1_MV^2)
            LossM1_sepLV <- loss_each_lv(model = model1_pls,
                                         PredErrors = PredErrorsM1_MV,
                                         compare_on_constructs = compare_on_constructs)
            LossM2 <- rowMeans(PredErrorsM2_MV^2)
            LossM2_sepLV <- loss_each_lv(model = model2_pls,
                                         PredErrors = PredErrorsM2_MV,
                                         compare_on_constructs = compare_on_constructs)
          }
  )

  # RMSE on indicator level
  switch (mod_comparison_type,
          m1_bench_m2_pls = {
            RMSE_M2 <- sqrt(colMeans(PredErrorsM2_MV^2))
            RMSE_IA <- sqrt(colMeans(PredErrors_IA^2))
            RMSE_LM <- sqrt(colMeans(PredErrors_LM^2))
            RMSE_mv <- cbind(RMSE_M2,
                             RMSE_LM,
                             RMSE_IA,
                             "Q^2_predict" = sapply(names(RMSE_M2), FUN = function(x){ 1-((RMSE_M2^2)[x])/((RMSE_IA^2)[x]) }))
          },
          m1_pls_m2_bench = {
            RMSE_M1 <- sqrt(colMeans(PredErrorsM1_MV^2))
            RMSE_IA <- sqrt(colMeans(PredErrors_IA^2))
            RMSE_LM <- sqrt(colMeans(PredErrors_LM^2))
            RMSE_mv <- cbind(RMSE_M1,
                             RMSE_LM,
                             RMSE_IA,
                             "Q^2_predict" = sapply(names(RMSE_M1), FUN = function(x){ 1-((RMSE_M1^2)[x])/((RMSE_IA^2)[x]) }))
          },
          m1_pls_m2_pls = {
            RMSE_M1 <- sqrt(colMeans(PredErrorsM1_MV^2))
            RMSE_M2 <- sqrt(colMeans(PredErrorsM2_MV^2))
            RMSE_M2 <- RMSE_M2[names(RMSE_M1)]
            RMSE_mv <- cbind(RMSE_M1,
                             RMSE_M2)
          }
  )
  # Output
  list(LossM1 = LossM1, LossM2 = LossM2, LossM1_sepLV = LossM1_sepLV, LossM2_sepLV = LossM2_sepLV, PredErrorsM1_MV = PredErrorsM1_MV, PredErrorsM2_MV = PredErrorsM2_MV,
       RMSE_mv = RMSE_mv)
}


# Helpers
check_and_classify_fun_input <- function(N,
                                         mv,
                                         cv_folds,
                                         model1,
                                         model2){
  # Cannot have both models as benchmark
  if (model1 == "benchmark" & model2 == "benchmark") {
    stop("At least one of the models must be a PLS path model")
  }
  # Classify which model is the benchmark model
  if (model1 == "IA" | model1 == "LM") {
    model1 <- "benchmark"
  }

  if (model2 == "IA" | model2 == "LM") {
    model2 <- "benchmark"
  }

  # Determine model comparison type
  mod_comparison_type <- c(m1_bench_m2_pls = FALSE,
                           m1_pls_m2_bench = FALSE,
                           m1_pls_m2_pls = FALSE)
  if (model1 == "benchmark" & model2 != "benchmark") {
    mod_comparison_type["m1_bench_m2_pls"] = TRUE
    }
  if (model1 != "benchmark" & model2 == "benchmark") {
    mod_comparison_type["m1_pls_m2_bench"] = TRUE
    }
  if (model1 != "benchmark" & model2 != "benchmark") {
    mod_comparison_type["m1_pls_m2_pls"] = TRUE
  }

  # Return
  names(mod_comparison_type[mod_comparison_type == TRUE])
}


est_pls_fun <- function(mv,
                        model){
  est_pls <- csem(mv,
                     model,
                     .disattenuate = FALSE)
# Ordering of variables
order_mv <- colnames(est_pls$Information$Data)
# Relation between MV and LV
MVLVrel <- est_pls$Information$Model$measurement
# Subset prediction errors for endogenous MV and reflective measurements
reflective_LV <- names(est_pls$Information$Model$construct_type[est_pls$Information$Model$construct_type == "Common factor"])
pred_relevant_LV <- intersect(est_pls$Information$Model$cons_endo,reflective_LV)
ReflEndogLV <- pred_relevant_LV
# Get the formative endogenous LV for later user message
FormativeLV <- names(est_pls$Information$Model$construct_type[est_pls$Information$Model$construct_type == "Composite"])
FormEndogLV <- intersect(est_pls$Information$Model$cons_endo, FormativeLV)

if (is.null(nrow(MVLVrel[ReflEndogLV,]))) {
  ReflEndogMV <- colnames(MVLVrel)[MVLVrel[ReflEndogLV,]!=0]
} else{
  pred_relevant_MV <- colnames(MVLVrel)[colSums(MVLVrel[ReflEndogLV,])==1]
  # Call them endogenous MV to be compatible with old code - must be change
  ReflEndogMV <- pred_relevant_MV
 }
# Exogenous variables model
PureExogLV <- est_pls$Information$Model$cons_exo
# Estimates
Loadings <- est_pls$Estimates$Loading_estimates
Weights <- est_pls$Estimates$Weight_estimates
Inner <- est_pls$Estimates$Path_estimates
# return value
list(est_pls = est_pls,
     order_mv = order_mv,
     MVLVrel = MVLVrel,
     ReflEndogLV = ReflEndogLV,
     FormEndogLV = FormEndogLV,
     ReflEndogMV = ReflEndogMV,
     PureExogLV = PureExogLV,
     Loadings = Loadings,
     Weights = Weights,
     Inner = Inner)
}

out_of_sample_pls <- function(model,
                              HoldoutMV,
                              FoldSize,
                              rescale_pred_errors,
                              TrainMV,
                              Train_std,
                              Train_mean){


  train_pls <- est_pls_fun(mv = TrainMV,
                           model = model)

  # Order holdout set according to PLS model
  HoldoutMV <- HoldoutMV[,train_pls$order_mv]
  HoldoutMV <- matrix(HoldoutMV,
                      ncol = length(train_pls$order_mv),
                      dimnames = list(NULL, train_pls$order_mv))
  # Calculate proxy of LV from indicators
  LVProxy<-HoldoutMV%*%t(train_pls$Weights)
  # Predict endog LV from exog LV
  PredictLV<-LVProxy%*%t(train_pls$Inner)
  # Predict MV connected to endogenous LV
  PredictMV<-PredictLV%*%train_pls$Loadings
  # Get pseudo prediction errors (i.e. also "prediction" errors on exogenous constructs)
  if (rescale_pred_errors == FALSE) {
    PseudoPredErrors<- HoldoutMV - PredictMV
  }
  if (rescale_pred_errors == TRUE) {
    # Rescale holdout
    HoldoutMV <- sapply(colnames(HoldoutMV),
                        FUN = function(x){HoldoutMV[,x]*Train_std[x] + Train_mean[x]},
                        USE.NAMES = FALSE)
    # Ensure that we have a matrix in the leave-one-out case
    HoldoutMV <- matrix(HoldoutMV,
                        ncol = length(train_pls$order_mv),
                        dimnames = list(NULL,train_pls$order_mv))
    # Rescale predictions
    PredictMV <- sapply(colnames(PredictMV),
                        FUN = function(x){PredictMV[,x]*Train_std[x] + Train_mean[x]},
                        USE.NAMES = FALSE)
    # Ensure that we have a matrix in the leave-one-out case
    PredictMV <- matrix(PredictMV,
                        ncol = length(train_pls$order_mv),
                        dimnames = list(NULL,train_pls$order_mv))

    PseudoPredErrors <- HoldoutMV - PredictMV
  }
  #PredErrors <- PseudoPredErrors[,train_pls$ReflEndogMV]
  PredErrors <- PseudoPredErrors[,train_pls$ReflEndogMV]
  PredErrors <- matrix(PredErrors,
                      ncol = length(train_pls$ReflEndogMV),
                      dimnames = list(NULL,train_pls$ReflEndogMV))


  # Output
  list(PredictMV = PredictMV,
       order_mv = train_pls$order_mv,
       order_endog_mv = train_pls$ReflEndogMV,
       PseudoPredErrors = PseudoPredErrors,
       PredErrors = PredErrors,
       train_pls = train_pls)
}
out_of_sample_bench <- function(out_samp_pls,
                                HoldoutMV,
                                FoldSize,
                                rescale_pred_errors,
                                TrainMV,
                                Train_std,
                                Train_mean){
  # Order the holdout set according to pls model
  HoldoutMV <- HoldoutMV[,out_samp_pls$order_mv]
  HoldoutMV <- matrix(HoldoutMV,
                      ncol = length(out_samp_pls$order_mv),
                      dimnames = list(NULL,out_samp_pls$order_mv))
  # Predictions: Indicator averages
  PredictMV_IA <- HoldoutMV*0
  lin_bench <- linear_bench_fun(out_samp_pls = out_samp_pls,
                   HoldoutMV = HoldoutMV,
                   TrainMV = TrainMV)


  if (rescale_pred_errors == FALSE) {
    PseudoPredErrors_IA <- HoldoutMV - PredictMV_IA
    PseudoPredErrors_LM <- HoldoutMV - lin_bench$PredictMV_LM
  }
  if (rescale_pred_errors == TRUE) {
    # Rescale holdout
    HoldoutMV <- sapply(colnames(HoldoutMV),
                        FUN = function(x){HoldoutMV[,x]*Train_std[x] + Train_mean[x]},
                        USE.NAMES = FALSE)
    HoldoutMV <- matrix(HoldoutMV,
                        ncol = length(out_samp_pls$order_mv),
                        dimnames = list(NULL, out_samp_pls$order_mv))
    # Rescale predictions, indicator average
    PredictMV_IA <- sapply(colnames(PredictMV_IA),
                        FUN = function(x){PredictMV_IA[,x]*Train_std[x] + Train_mean[x]},
                        USE.NAMES = FALSE)
    PredictMV_IA <- matrix(PredictMV_IA,
                        ncol = length(out_samp_pls$order_mv),
                        dimnames = list(NULL, out_samp_pls$order_mv))
    # Pseudo prediction errors, indicator average
    PseudoPredErrors_IA <- HoldoutMV - PredictMV_IA
    # Rescale predictions, linear model
    lin_bench$PredictMV_LM <- sapply(colnames(lin_bench$PredictMV_LM),
                           FUN = function(x){lin_bench$PredictMV_LM[,x]*Train_std[x] + Train_mean[x]},
                           USE.NAMES = FALSE)
    lin_bench$PredictMV_LM <- matrix(lin_bench$PredictMV_LM,
                           ncol = length(out_samp_pls$order_mv),
                           dimnames = list(NULL, out_samp_pls$order_mv))
    # Pseudo prediction errors, linear model
    PseudoPredErrors_LM <- HoldoutMV - lin_bench$PredictMV_LM
  }
  PredErrors_IA <- PseudoPredErrors_IA[,out_samp_pls$order_endog_mv]
  PredErrors_LM <- PseudoPredErrors_LM[,out_samp_pls$order_endog_mv]
  # Output
  list(PredictMV_IA = PredictMV_IA,
       PseudoPredErrors_IA = PseudoPredErrors_IA,
       PredErrors_IA = PredErrors_IA,
       PredictMV_LM = lin_bench$PredictMV_LM,
       PseudoPredErrors_LM = PseudoPredErrors_LM,
       PredErrors_LM = PredErrors_LM,
       order_mv = out_samp_pls$order_mv)
}

linear_bench_fun <- function(out_samp_pls, HoldoutMV, TrainMV){
  # Identify exogenous MV for each endogenous MV
  # 1. Identify MV that is exogenous to each LV
  MV_exog_to_LV <- out_samp_pls$train_pls$est_pls$Information$Model$structural %*%
    out_samp_pls$train_pls$est_pls$Information$Model$measurement
  # 2. Identify exogenous MV for each endogenous MV
  MV_exog_to_MV <- t(t(MV_exog_to_LV) %*% out_samp_pls$train_pls$est_pls$Information$Model$measurement)
  # Predictions: Linear benchmark
  PredictMV_LM <- HoldoutMV*0
  # Empty list to hold results
  train_lm <- vector(mode = "list", length = length(out_samp_pls$order_endog_mv))
  names(train_lm) <- out_samp_pls$order_endog_mv
  for (k in out_samp_pls$order_endog_mv) {
    train_lm[[k]] <- lm(data = as.data.frame(TrainMV),
                        formula = as.formula(paste(c(k,
                                                     "~" ,
                                                     paste(colnames(MV_exog_to_MV)[MV_exog_to_MV[k,]==1], collapse = "+")
                        ),
                        collapse=" ")))
    PredictMV_LM[,k] <- stats::predict(train_lm[[k]], newdata = as.data.frame(HoldoutMV))
  }
  return(list(MV_exog_to_LV = MV_exog_to_LV,
              MV_exog_to_MV = MV_exog_to_MV,
              PredictMV_LM = PredictMV_LM,
              train_lm = train_lm))
}

loss_each_lv <- function(model,
                         PredErrors,
                         compare_on_constructs){
  Loss_sepLV <- vector(mode = "list", length = length(model$ReflEndogLV))
  names(Loss_sepLV) <- model$ReflEndogLV
  for (j in model$ReflEndogLV) {
    if (length(colnames(model$MVLVrel)[model$MVLVrel[j,]==1])>1) {
    Loss_sepLV[[j]]<-rowMeans(PredErrors[,colnames(model$MVLVrel)[model$MVLVrel[j,]==1]]^2)
    }else{
      Loss_sepLV[[j]]<-PredErrors[,colnames(model$MVLVrel)[model$MVLVrel[j,]==1]]^2
    }
  }

  if (is.character(compare_on_constructs)) {
    comp_LV <- paste(compare_on_constructs, collapse = "_")
    # If we compare on more than one construct, we will have multiple indicators (use rowMeans)
    if (length(compare_on_constructs) > 1) {
      ## Losses for set of constructs
      Loss_sepLV[[comp_LV]] <- rowMeans(PredErrors[,colnames(model$MVLVrel)[(colSums(model$MVLVrel[compare_on_constructs,]==1)==1)]]^2)
    }
    # If we compare on a single construct, then the loss for this construct is already calculated in Loss_sepLV
    if (length(compare_on_constructs) == 1) {
      # If we compare on a single construct, then the loss for this construct is already calculated in Loss_sepLV
      # TODO: Should we give a warning here?
    }

  }
  Loss_sepLV
}
