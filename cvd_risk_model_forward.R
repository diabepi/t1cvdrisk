)###############################################################################
#  Created on 1 Jun 2020
#  ------------------------------------------------------------
#  Copyright (c) 2020 Diabetes Epidemiology Group.
#  All Right Reserved.
#  ------------------------------------------------------------
#  Author:  smcgurnaghan

# rm(list=ls()) 

svn_author  = "$Author: stuart $"
svn_version = "$Rev: 9323 $"
svn_date    = "$Date: 2020-11-25 18:00:52 +0000 (Wed, 25 Nov 2020) $"
svn_id      = "$Id: cvd_risk_model_forward.R 9323 2020-11-25 18:00:52Z stuart $"
svn_URL     = "$URL: https://diabepi.igmm.ed.ac.uk/svn/sdrnepi/trunk/projects/national_dataset/src/t1cvdrisk/cvd_risk_model_forward.R $"

#' @description 
#' Build a model using selection to select the best predictors of the
#' defined outcome. Perform cross-fold validation of the resulting model
#' @param 
#' opts options object containing all study options
#' 
#' @param
#' surv.obj the survival dataset where we only use the exit slice 
#' 
#' 
cvdRiskModelForward <- function(opts, surv.obj) {
  
  library(wevid)
  
  outcome <- opts$options()$model.outcome
  folds.n <- opts$options()$model.folds
  
  # Minimum model of age, sex and diabetes duration
  basecols.rgx <- opts$options()$model.min
  
  base.cols <- names(
      surv.obj$retcols.regex(
          basecols.rgx,
          baseline = TRUE
      )
  )
  
  # Full model of all predictors / independent variables
  fullcols.rgx <- opts$options()$model.candidates
  
  full.cols <- names(
      surv.obj$retcols.regex(
          paste0(
              fullcols.rgx
          ),
          baseline = TRUE
      )
  )
  
  # If using MI, substitute factors with 'all chain' weighted design matrices
  if(opts$options()$impute.baseline) {
    full.cols <- gsub("simd", "simd.mm", full.cols)
    full.cols <- gsub("smoker.evr.max", "smoker.mm", full.cols)
    full.cols <- gsub("albuminuric.grade.max", "albuminuric.grade.mm", full.cols)
    full.cols <- gsub("retgrade", "retgrade.mm", full.cols)
    full.cols <- gsub("ethnic.sub", "ethnic.sub.mm", full.cols)
  }
  
  # Save models info for documentation
  saveRDS(full.cols, dataPath("full.cols.forward.RData"))
  saveRDS(base.cols, dataPath("base.cols.forward.RData"))
  saveRDS(with(xx <- data.table(full.cols), 
          print(varsToStrings(xx,"full.cols"))), dataPath("full.cols.forward.txt.RData"))
  
  # Finalise model dataset
  dataset <- na.omit(surv.obj$retcols.regex(
          paste0(paste(full.cols,collapse="|"), 
              "|slice|entry.age|^age.|sex|duration|slice.fu|year|", outcome)))
  
  # Retrieve imputed average model matrix data for factors
  if(opts$options()$impute.baseline) {
    simd.mm <- as.model.matrix(na.omit(surv.obj$retcols.regex("simd")), 
        "simd", c('SIMD2','SIMD3','SIMD4','SIMD5'))[,-1]
    albuminuric.grade.mm <- 
        as.model.matrix(na.omit(surv.obj$retcols.regex("albuminuric.grade.max")), 
            "albuminuric.grade.max", c('Micro','Macro'))[,-1]
    retgrade.mm <- as.model.matrix(na.omit(surv.obj$retcols.regex("retgrade")), 
        "retgrade", c('Non Referable', 'Referable / Eye Clinic'))[,-1]
    smoker.mm <- as.model.matrix(na.omit(surv.obj$retcols.regex("smoker.evr.max")), 
        "smoker.evr.max", c('Smoker'))[,-1]
  }
  
  # Narrow the model columns for adjustment - only continuous values
  adj.cols <- unique(c(full.cols[grepl("med|avg|crt", full.cols)],
      base.cols[grepl("med|avg|crt", base.cols)]))
  
  
  # Add squared terms to continuous data and squared logs for skewed data
  if(opts$options()$check.skewness) {
    
    var.skewed = ""
    var.remove = as.character()
    
    for (modvar in unique(c(adj.cols))) {
      
      if(is.numeric(surv.obj$surv.table[[modvar]]) & 
          length(unique((surv.obj$surv.table[[modvar]]))) > 20 & 
          modvar != "slice.fu") 
      {
       
        if(abs(e1071::skewness(surv.obj$surv.table[[modvar]]))>1) {
          R.utils::printf("%s is skewed in set \n", modvar, var.skewed)
          var.skewed= paste(var.skewed, paste0(" + log(", modvar, ") "))
          # Remove original vars as they have been adjusted for skewness 
          var.remove = c(var.remove, modvar)
        } 
      }
    }  
   var.skewed= gsub("\\+ * \\+", "\\+", var.skewed)
  } else {
   var.skewed= ""
  }
  
  var.remove <- c(var.remove, "egfr.med")
  adj.cols <- setdiff(adj.cols, var.remove)
  full.cols <- setdiff(full.cols, var.remove)
  
  # Add interactions with age and sex
  if(opts$options()$interactions.sex) {
    # Additional manually added post hoc
    adj.cols.s <- c(adj.cols, "log(egfr.med)", "log(bmi.med)", "log(tchol.hdl.med)")
    intact.sex <- paste0(paste0(adj.cols.s[!grepl("sex|age.crt", adj.cols)], " * sex"), " + ")
  } else {
    intact.sex = ""
  }
  
  if(opts$options()$interactions.age) {
    # Additional manually added post hoc
    adj.cols.a <- c(adj.cols, "log(egfr.med)", "log(bmi.med)", "log(tchol.hdl.med)")
    intact.age <- paste0(paste0(adj.cols.a[!grepl("age.crt", adj.cols)], " * age.crt"), " + ")
  } else {
    intact.age = ""
  }
  
  # Base model with offset term
  base.model <- glm(
      formula(paste0(outcome, " ~  offset(log(slice.fu)) + ",  
              paste(base.cols, collapse="+"))), 
      family="poisson",
      data = na.omit(dataset)
  )
  
  # Full model with offset term, interactions, squared terms and logs for skewed data
  formula.in <- gsub("\\+ * \\+","+", 
      paste0(outcome, " ~  offset(log(slice.fu)) + ",  
          paste(intact.age, sep="", collapse=""),
          paste(intact.sex, sep="", collapse=""),
          var.skewed, " + eval(age.crt^3) + eval(age.crt^2) + log(egfr.med) + ",
          paste(paste(unique(c(base.cols, full.cols)), collapse="+"))))
  
  full.model <- glm(
      formula(formula.in), 
      family="poisson",
      data = na.omit(dataset)
  )
  
  # Run selection
  set.seed(opts$options()$seed)
  forward.model <-
      createOrRetrieve(
          "forward.model.Rdata",
          creator = MASS::stepAIC,
          creator.args = list(
              object = base.model,
              scope = list(lower=base.model, upper=full.model),
              direction = "forward",
              trace = 1,
              k = 3,
              method = "approximate"
          ), skip.global = TRUE 
      )
  
  
  
  
  
  # Crude C-statistic as sanity.check
  pROC::auc(dataset[[outcome]], 
      predict(forward.model, type="response") / dataset$slice.fu)
  
  
  
  
  # Crude calibration plot as sanity.check
  rDiabSource("functions.R", topic = "common", 
      project = "national_dataset", local=F) 
  library(Hmisc)
  library(gmodels)
  pdf(file = repPath("forward_model_sanity.pdf"))
  forward.crude <- calibration.plot.cex(dataset[,.(
              plotID=.I,
              Observed = get(outcome),
              Predicted1 = predict(forward.model, type="response")/slice.fu)],
      n.years = 1,
      which.model = 1,
      ci.lines = TRUE,
      na.rm = TRUE,
      alpha = 0.05,
      N.bins = 10,
      xlab = "Predicted Risk of CVD (deciles)",
      ylab = "Observed risk as proportion of patients",
      model.names = NULL,
      main = paste0("All crude (",
          nrow(dataset)," person/yrs)"))
  dev.off()
  

  
  dataset.incprior <- rbind(surv.table.priors, surv.obj$surv.table)
  
  
  # Bring in the model
  model.in <- readRDS(dataPath(paste0("forward.model.Rdata")))$data.obj
  # Extract the environment
  model.env <- environment(model.in$formula)
  # Update the model env data
  model.env$simd.mm <- 
      as.model.matrix(dataset.incprior[,.(simd)], 
          "simd", c('SIMD2','SIMD3','SIMD4','SIMD5'))[,-1]
  model.env$albuminuric.grade.mm <- 
      as.model.matrix(dataset.incprior[,.(albuminuric.grade.max)], 
          "albuminuric.grade.max", c('Micro','Macro'))[,-1]
  model.env$retgrade.mm <- 
      as.model.matrix(dataset.incprior[,.(retgrade)],
          "retgrade", c('Non Referable', 'Referable / Eye Clinic'))[,-1]
  model.env$smoker.mm <- 
      as.model.matrix(dataset.incprior[,.(smoker.evr.max)], 
          "smoker.evr.max", c('Smoker'))[,-1]
  

  
  # Crude C-statistic as sanity.check

pROC::auc(dataset.incprior[[outcome]], 
    predict(model.in, newdata = dataset.incprior, type="response") / dataset.incprior$slice.fu)


  
  # Crude calibration plot as sanity.check + priors
  rDiabSource("functions.R", topic = "common", 
      project = "national_dataset", local=F) 
  library(Hmisc)
  library(gmodels)
  pdf(file = repPath("forward_model_sanity_incprior.pdf"))
  forward.crude <- calibration.plot.cex(dataset.incprior[,.(
              plotID=.I,
              Observed = get(outcome),
              Predicted1 = predict(forward.model, 
                  newdata=dataset.incprior,
                  type="response")/dataset.incprior$slice.fu)],
      n.years = 1,
      which.model = 1,
      ci.lines = TRUE,
      na.rm = TRUE,
      alpha = 0.05,
      N.bins = 10,
      xlab = "Predicted Risk of CVD (deciles)",
      ylab = "Observed risk as proportion of patients",
      model.names = NULL,
      main = paste0("All crude (",
          nrow(dataset.incprior)," person/yrs)"))
  dev.off()
  
  
  
  
  
  forward.model.table <- 
      data.table::data.table(jtools::summ(forward.model, 
              confint=TRUE, exp=TRUE, )$coeftable, keep.rownames = TRUE)
  
  data.table::data.table(summary(forward.model)$coefficients, keep.rownames = TRUE)
  data.table::setnames(forward.model.table, "rn", "Predictor")
  data.table::setnames(forward.model.table,  "p", "P-Value")
  data.table::setnames(forward.model.table,  "exp(Est.)", "IRR")
  varsToStrings(forward.model.table, "Predictor")
  forward.model.table <- forward.model.table[,-5]
  saveRDS(forward.model.table, dataPath("forward.model.table.Rdata"))
  
  # Remove cohort information for model export
  forward.model.redact <- data.table::copy(forward.model)
  forward.model.redact[c("residuals", "weights", "fitted.values","data","prior.weights",
          "na.action","effects","y","offset","predictors",
          "linear.predictors","model")] <- NULL
  save(forward.model.redact, file=dataPath("scottish.forward.model.RData"))
  
  # Prediction performance
  dataset <- surv.obj$retcols.regex(
          paste0(paste(full.cols,collapse="|"), 
              "|serialno|age.crt|entry|year|sex|*.med|*.avg|duration|",
              "retgrade|simd|albuminuric|smoker.evr.max|slice.fu|", outcome)) 
  stopifnot(nrow(na.omit(dataset)) == nrow(surv.obj$surv.table))
  
  
  set.seed(opts$options()$seed)
  fold.indices <- caret::createFolds(dataset[[outcome]], k = folds.n)
  fold.results.test <- fold.results.train <- full.densities <- list()
  fold.auc <- vector()
  
  #' Overload the step program to use the createOrRetrieve
  #' 
  #' 
  stepw  <- function(traindata = traindata, ...) { 
    traindata = traindata
    return(step(...));
  }
  
  #' Perform cross fold regression - parallel capable function
  #' with backup loading via createOrRetrieve
  #'
  stepFold <- function(i, traindata) {
    
    fold <- names(fold.indices[i])
    
    # Retrieve imputed average model matrix data for factors
    simd.mm <- as.model.matrix(traindata[,.(simd)], 
        "simd", c('SIMD2','SIMD3','SIMD4','SIMD5'))[,-1]
    albuminuric.grade.mm <- as.model.matrix(traindata[,.(albuminuric.grade.max)], 
        "albuminuric.grade.max", c('Micro','Macro'))[,-1]
    retgrade.mm <- as.model.matrix(traindata[,.(retgrade)], 
        "retgrade", c('Non Referable', 'Referable / Eye Clinic'))[,-1]
    smoker.mm <- as.model.matrix(traindata[,.(smoker.evr.max)], 
        "smoker.evr.max", c('Smoker'))[,-1]
    
    
    # Base model with offset term
    base.model <- glm(
        formula(paste0(outcome, " ~  offset(log(slice.fu)) + ",  
                paste(base.cols, collapse="+"))), 
        family="poisson",
        data = traindata
    )
    saveRDS(base.model, dataPath("base.model.Rdata"))
    
    
    # Full model with offset term
    full.model <- glm(
        formula(paste0(outcome, " ~  offset(log(slice.fu)) + ",  
                paste(intact.age, sep="", collapse=""),
                paste(intact.sex, sep="", collapse=""),
                var.skewed, " + ",
                paste(paste(c(base.cols, full.cols), collapse="+")))), 
        family="poisson",
        data = traindata
    )
    
    set.seed(opts$options()$seed)
    fold.model.train <-
        createOrRetrieve(
            sprintf("forward.%s.Rdata", fold),
            creator = stepw,
            creator.args = list(
                object = base.model,
                scope = list(lower=base.model, upper=full.model),
                trace = -1,
                k = 3,
                direction = "forward",
                method = "approximate",
                traindata = traindata
            ), skip.global = TRUE
        )
    return(fold.model.train)
  }   
  
  # Scope is global for the results
  fold.model.train <- list()
  fold.model.test <- list()
  
  # Parallelised
  library(doMC)
  registerDoMC(cores=10)
  
  fold.model.train <- foreach(i = 1:folds.n) %dopar% {
    
    stepFold(i, dataset[!fold.indices[[i]]])
    
  }
  
  # Predictive performance - testing set folds
  foldMod <- wevid.fold <- list()
  for (nf in 1:folds.n) { 
    
    nf <- sprintf("%02d",nf)  
    
    # Bring in the model
    model.in <- readRDS(dataPath(paste0("forward.Fold", nf,".Rdata")))$data.obj
    # Extract the environment
    model.env <- environment(model.in$formula)
    # Update the model env data
    model.env$simd.mm <- 
        as.model.matrix(dataset[fold.indices[[paste0("Fold",nf)]]][,.(simd)], 
            "simd", c('SIMD2','SIMD3','SIMD4','SIMD5'))[,-1]
    model.env$albuminuric.grade.mm <- 
        as.model.matrix(dataset[fold.indices[[paste0("Fold",nf)]]][,.(albuminuric.grade.max)], 
            "albuminuric.grade.max", c('Micro','Macro'))[,-1]
    model.env$retgrade.mm <- 
        as.model.matrix(dataset[fold.indices[[paste0("Fold",nf)]]][,.(retgrade)],
            "retgrade", c('Non Referable', 'Referable / Eye Clinic'))[,-1]
    model.env$smoker.mm <- 
        as.model.matrix(dataset[fold.indices[[paste0("Fold",nf)]]][,.(smoker.evr.max)], 
            "smoker.evr.max", c('Smoker'))[,-1]

    
    # Predict probability in testing folds (loads model for fold)
    # The predict function automatically calculates an adjusted probability
    # for intervals of length < 1.
    prob.pred.raw <- predict.glm( model.in, 
        newdata=dataset[fold.indices[[paste0("Fold",nf)]]], type="response")
    
    # With Poisson, we adjust probability back to equivalent intervals for calculating
    # the AUC / C-statistic
    prob.pred <- prob.pred.raw / dataset[fold.indices[[paste0("Fold",nf)]]]$slice.fu
    
    
    # Training fold prior probability of an event
    prior.prob.calc <- mean(dataset[!fold.indices[[paste0("Fold",nf)]]][[outcome]]) 
    
    dataset[fold.indices[[paste0("Fold",nf)]], fold := nf]
    # Training population prior probability 
    dataset[fold.indices[[paste0("Fold",nf)]], prior.prob := prior.prob.calc]
    # Test population predicted probability 
    dataset[fold.indices[[paste0("Fold",nf)]], pred.prob := prob.pred]
    # Store the unadjusted probability
    dataset[fold.indices[[paste0("Fold",nf)]], pred.prob.raw := prob.pred.raw]  
    
    # Check for log hazard outwith the floating point of the system
    stopifnot(nrow(dataset[fold.indices[[paste0("Fold",nf)]]][pred.prob<=0 | pred.prob >= 1])<1)
   
    # Store the result for each fold in the event we need to analyse
    wevid.fold[[nf]] <-  
        Wdensities(
            y = dataset[fold.indices[[paste0("Fold",nf)]]][[outcome]], 
            posterior.p = dataset[fold.indices[[paste0("Fold",nf)]]]$pred.prob, 
            prior.p = dataset[fold.indices[[paste0("Fold",nf)]]]$prior.prob, 
            recalibrate = FALSE
        )
  }
  
  # Weights of evidence across all testing folds (100% of the data)
  wevid.result.forward <-   Wdensities(
      y = dataset[[outcome]], 
      posterior.p = dataset$pred.prob, 
      prior.p = dataset$prior.prob,  
      recalibrate = FALSE
  )
  
  pdf(file=paste0(repPath("wevid_forward_full_densities.pdf")))
  plotWdists(wevid.result.forward)
  dev.off()
  pdf(file=paste0(repPath("wevid_forward_full_roc.pdf")))
  plotroc(wevid.result.forward)
  dev.off()
  
  pROC::roc(dataset[[outcome]], dataset$pred.prob)
  pdf(file=paste0(repPath("roc_forward_plot.pdf")))
  pROC::plot.roc(dataset[[outcome]], dataset$pred.prob, legacy.axes = TRUE)
  dev.off()
  
  
  wevid.result.forward.tab <- 
      data.table(t(summary(wevid.result.forward)), keep.rownames = TRUE)
  names(wevid.result.forward.tab) <- c("Parameter", "Full Model")
  
  
  # Calculate the increment in prediction from the base
  training.base.model.folds <- list()
  for(i in 1:folds.n) {
    training.base.model.folds[[i]] <- glm(
        formula(paste0(outcome, " ~ ",  
                paste(base.cols, collapse="+"))), 
        family="binomial",
        data = na.omit(dataset[-fold.indices[[i]]])
    )
  }
  
  
  # Predictive performance - testing folds 
  foldModbase <- wevid.foldbase <- list()
  for (nf in 1:folds.n) { 

    # Calculate predictions across all folds for all testees.
    foldpredBase <- predict(training.base.model.folds[[nf]], 
        newdata=dataset[fold.indices[[nf]]],  type="response")
    dataset[fold.indices[[nf]], pred.base := foldpredBase]
    # Training fold prior probability of an event 
    prior.prob.calc <- mean(dataset[!fold.indices[[nf]]][[outcome]])   
    # Training population prior probability 
    dataset[fold.indices[[nf]], prior.base := prior.prob.calc]
  }
  
  wevid.forward.base <- 
      Wdensities(
          y = dataset[[outcome]], 
          posterior.p = dataset$pred.base, 
          prior.p = dataset$prior.base, 
          recalibrate = FALSE
      )
  
  pdf(file=paste0(repPath("wevid_forward_base_densities.pdf")))
  plotWdists(wevid.forward.base)
  dev.off()
  pdf(file=paste0(repPath("wevid_forward_base_roc.pdf")))
  plotroc(wevid.forward.base)
  dev.off()
  
  wevid.forward.base.tab <- 
      data.table(t(summary(wevid.forward.base)), keep.rownames = TRUE)
  names(wevid.forward.base.tab) <- c("Parameter", "Base Model")
  
  wevid.forward.result.tab <- 
      cbind(wevid.forward.base.tab, wevid.result.forward.tab[,2])
  saveRDS(wevid.forward.result.tab, dataPath("wevid.forward.results.tab.Rdata"))
  
  saveRDS(pROC::ci.auc(dataset[[outcome]], dataset$pred.base), 
      dataPath("confint.forward.auroc.base.Rdata"))
  auroc.forward.ci.base <- readRDS(dataPath("confint.forward.auroc.base.Rdata"))
  
  saveRDS(pROC::ci.auc(dataset[[outcome]], dataset$pred.prob), 
      dataPath("confint.forward.auroc.full.Rdata"))
  auroc.forward.ci.full <- readRDS(dataPath("confint.forward.auroc.full.Rdata"))
  

  forward.model.table <- readRDS(dataPath("forward.model.table.Rdata"))
  wevid.result.forward.tab <- readRDS(dataPath("wevid.forward.results.tab.Rdata"))
  
  saveRDS(dataset, dataPath("forward_model_predictions.Rdata"))
  
  # Add age bands
  dataset[, roc.age.band := cut(age.crt,breaks=c(20,40,60,+Inf))]
  levels(dataset$roc.age.band) <- c("20-40","40-60","60+")
  
  wevid.result.1yr.female <- with(dataset[sex=="Female"],
      Wdensities(
          cvd.any.bin,
          pred.prob,
          prior.prob
      )
  )
  
  wevid.result.1yr.male <- with(dataset[sex=="Male"],
      Wdensities(
          cvd.any.bin,
          pred.prob,
          prior.prob
      )
  )
  
  
  wevid.result.1yr.20to40 <- with(dataset[roc.age.band=='20-40'],
      Wdensities(
          cvd.any.bin,
          pred.prob,
          prior.prob
      )
  )
  
  wevid.result.1yr.40to60 <- with(dataset[roc.age.band=='40-60'],
      Wdensities(
          cvd.any.bin,
          pred.prob,
          prior.prob
      )
  )
  
  wevid.result.1yr.60P <- with(dataset[roc.age.band=='60+'],
      Wdensities(
          cvd.any.bin,
          pred.prob,
          prior.prob
      )
  )
  
  wevid.result.1yrb.female <- with(dataset[sex=="Female"],
      Wdensities(
          cvd.any.bin,
          pred.base,
          prior.base
      )
  )
  
  wevid.result.1yrb.male <- with(dataset[sex=="Male"],
      Wdensities(
          cvd.any.bin,
          pred.base,
          prior.base
      )
  )
  
  wevid.result.1yrb.20to40 <- with(dataset[roc.age.band=='20-40'],
      Wdensities(
          cvd.any.bin,
          pred.base,
          prior.base
      )
  )
  
  wevid.result.1yrb.40to60 <- with(dataset[roc.age.band=='40-60'],
      Wdensities(
          cvd.any.bin,
          pred.base,
          prior.base
      )
  )
  
  wevid.result.1yrb.60P <- with(dataset[roc.age.band=='60+'],
      Wdensities(
          cvd.any.bin,
          pred.base,
          prior.base
      )
  )
  
  

  # Wevid table results	
  all.wevid.full <- data.table(rbind(
          cbind(Group="All", summary(wevid.result.forward)),
          cbind(Group="Male", summary(wevid.result.1yr.male)),
          cbind(Group="Female", summary(wevid.result.1yr.female)),
          cbind(Group="20-40", summary(wevid.result.1yr.20to40)),
          cbind(Group="40-60", summary(wevid.result.1yr.40to60)),
          cbind(Group="60+", summary(wevid.result.1yr.60P))
      ))
  all.wevid.full[,c(1,2,4,5,6,7)]	
  setnames(all.wevid.full, 
      'log-likelihood after recalibration (nats)', 'log-ll after recal (nats)')
  
  
  all.wevid.base <- data.table(rbind(
          cbind(Group="All", summary(wevid.forward.base)),
          cbind(Group="Male", summary(wevid.result.1yrb.male)),
          cbind(Group="Female", summary(wevid.result.1yrb.female)),
          cbind(Group="20-40", summary(wevid.result.1yrb.20to40)),
          cbind(Group="40-60", summary(wevid.result.1yrb.40to60)),
          cbind(Group="60+", summary(wevid.result.1yrb.60P))
      ))
  setnames(all.wevid.base, 
      'log-likelihood after recalibration (nats)', 'log-ll after recal (nats)')
  all.wevid.base[,c(1,4,6,7,8)]	
  
# Combined
  all.wevid.combined <- merge(all.wevid.full[,c(1,3,5,7,8)], 
      all.wevid.base[,c(1,3,5,7,8)], by="Group")
  all.wevid.combined <- all.wevid.combined[, 
      .(Group, `Crude C-statistic` = `Crude C-statistic.x` - 
              `Crude C-statistic.y`,
          `Crude Λ (bits)` = `Crude Λ (bits).x` - `Crude Λ (bits).y`,
          `Test log-likelihood (nats)` = `Test log-likelihood (nats).x` - 
              `Test log-likelihood (nats).y`,
          `log-likelihood after recal (nats)` = 
              `log-ll after recal (nats).x` - 
              `log-ll after recal (nats).y`)]
  
  
# Cut down side by side
  
  all.wevid.side <- merge(all.wevid.full[,c(1,3,5)], 
      all.wevid.base[,c(1,3,5)], all.x=T, by="Group")
  all.wevid.side <- all.wevid.side[,c(1,2,4,3,5)]
  all.wevid.side
  setnames(all.wevid.side, "Crude C-statistic.x", 
      "C-statistic (Final)")
  setnames(all.wevid.side, "Crude C-statistic.y", 
      "C-statistic (Base)")
  setnames(all.wevid.side, "Crude Λ (bits).x", 
      "Λ (Final)")
  setnames(all.wevid.side, "Crude Λ (bits).y", 
      "Λ (Base)")
# Format and order columns for paper
  all.wevid.side <- cbind(all.wevid.side, abs(all.wevid.base[,7]-all.wevid.full[,7]))
  fmat <- function(x) {return(format(round(x,2),snmall=2))}
  all.wevid.side <- all.wevid.side[, c(Group=list(Group),lapply(.SD, fmat)), .SDcols=c(2:6)]
  all.wevid.side <- all.wevid.side[,c(1,3,2,5,4,6)]
  setnames(all.wevid.side,"Test log-likelihood (nats)", "log-likelihood increment (natural log units)")
    
  saveRDS(all.wevid.side, dataPath("all.wevid.side.Rdata"))

  
  calibration <- function(dataset) {
    
    rDiabSource("functions.R", topic = "common", 
         project = "national_dataset", local=F) 
    library(Hmisc)
    library(gmodels)
        
    pdf(file=repPath("calibration_forward.pdf"))
    calib <- calibration.plot.cex(dataset[,.(plotID=.I,
                Observed=dataset$cvd.any.bin,Predicted1=pred.prob)],
        n.years=1,
        which.model=1,
        ci.lines=TRUE,
        na.rm=TRUE,
        alpha=0.05,
        N.bins=10,
        xlab="Predicted risk of incident CVD by decile of risk",
        ylab="Observed risk as proportion of patients",
        model.names=NULL
        )
    dev.off()
    
    pdf(file=repPath("failure.10year2040.pdf"))
    failure.10year2040 <- calibration.plot.cex(dataset[roc.age.band=='20-40',.(plotID=.I,
                Observed=cvd.any.bin,Predicted1=pred.prob)],
        n.years=1,
        which.model=1,
        ci.lines=TRUE,
        na.rm=TRUE,
        alpha=0.05,
        N.bins=10,
        xlab="Predicted Risk of CVD (deciles)",
        ylab="Observed risk as proportion of patients",
        model.names=NULL,
        main=paste0("Age at entry 20-40 (",
            nrow(dataset[roc.age.band=='20-40'])," person/yrs)"))
    dev.off()
    

    pdf(file=repPath("failure.10year4060.pdf"))
    failure.10year4060 <- calibration.plot.cex(dataset[roc.age.band=='40-60',.(plotID=.I,
                Observed=cvd.any.bin,Predicted1=pred.prob)],
        n.years=1,
        which.model=1,
        ci.lines=TRUE,
        na.rm=TRUE,
        alpha=0.05,
        N.bins=10,
        xlab="Predicted Risk of CVD (deciles)",
        ylab="Observed risk as proportion of patients",
        model.names=NULL,
        main=paste0("Age at entry 40-60 (",
            nrow(dataset[roc.age.band=='40-60'])," person/yrs)"))
    dev.off()
    
    
    pdf(file=repPath("failure.10year60P.pdf"))
    failure.10year60P <- calibration.plot.cex(dataset[roc.age.band=='60+',.(plotID=.I,
                Observed=cvd.any.bin,Predicted1=pred.prob)],
        n.years=1,
        which.model=1,
        ci.lines=TRUE,
        na.rm=TRUE,
        alpha=0.05,
        N.bins=10,
        xlab="Predicted Risk of CVD (deciles)",
        ylab="Observed risk as proportion of patients",
        model.names=NULL,
        main=paste0("Age at entry 60+ (",
            nrow(dataset[roc.age.band=='60+'])," person/yrs)"))
    dev.off()
    

    pdf(file=repPath("failure.10yearfem.pdf"))
    failure.10year.fem <- calibration.plot.cex(dataset[sex=="Female",.(plotID=.I,
                Observed=cvd.any.bin,Predicted1=pred.prob)],
        n.years=1,
        which.model=1,
        ci.lines=TRUE,
        na.rm=TRUE,
        alpha=0.05,
        N.bins=10,
        xlab="Predicted Risk of CVD (deciles)",
        ylab="Observed risk as proportion of patients",
        model.names=NULL,
        main=paste0("Female (",
            nrow(dataset[sex==1])," person/yrs)"))
    dev.off()
        
    
    pdf(file=repPath("failure.10yearmal.pdf"))
    failure.10year.mal <- calibration.plot.cex(dataset[sex=="Male",.(plotID=.I,
                Observed=cvd.any.bin,Predicted1=pred.prob)],
        n.years=1,
        which.model=1,
        ci.lines=TRUE,
        na.rm=TRUE,
        alpha=0.05,
        N.bins=10,
        xlab="Predicted Risk of CVD (deciles)",
        ylab="Observed risk as proportion of patients",
        model.names=NULL,
        main=paste0("Male (",
            nrow(dataset[sex==0])," person/yrs)"))
    dev.off()
    
    hoslem.forward <- with(dataset[,.(plotID=.I,
                Observed=cvd.any.bin,Predicted1=pred.prob)], 
        ResourceSelection::hoslem.test(Observed, Predicted1))
    saveRDS(hoslem.forward, dataPath("hoslem.forward.RData"))
  }
  
  # Calibration usinh the test dataset results
  calibration(dataset)
  
  # reinstate surv obj
  surv.obj <- readRDS(dataPath("surv.obj.Rdata"))$data.obj
  
}


