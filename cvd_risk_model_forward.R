##############################################################################
#  Created on 1 Jun 2020
#  ------------------------------------------------------------
#  Copyright (c) 2020 Diabetes Epidemiology Group.
#  All Right Reserved.
#  ------------------------------------------------------------
#  Author:  smcgurnaghan

# rm(list=ls()) 

svn_author  = "$Author$"
svn_version = "$Rev$"
svn_date    = "$Date$"
svn_id      = "$Id$"
svn_URL     = "$URL$"

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
cvdRiskModelForward <- function(
		opts, 
		surv.obj, 
		interval.model = opts$options()$gap, 
		interval.short = 14
) {
	
	library(wevid)
	library(MASS)
	library(ggplot2)
	
	
	
	
	test.calibration.table <- function(table.bins) {
		obs <- table.bins[1, ]
		exp <- table.bins[2, ]
		df <- ncol(table.bins)
		teststat.chisq <- sum((obs - exp) ^2 / exp)
		pvalue <- 1 - pchisq(q=teststat.chisq, df=df)
		return(c(teststat.chisq=teststat.chisq, df=df,pvalue=pvalue)) 
	}
	
	
	plot.calibration <- function(observed, expected, title.in = "") {
		
		riskdecile.group <- cut(expected, quantile(x=expected, probs=seq(0, 1, by=0.1)))
		expected.riskdecile <- tapply(expected, riskdecile.group, sum)
		observed.riskdecile <- tapply(observed, riskdecile.group, sum)
		calib <- data.frame(x=expected.riskdecile, y=observed.riskdecile)
		max.xy <- max(c(calib$x, calib$y))
		
		p <- ggplot(data=calib, aes(x=x, y=y)) + geom_point(size=2) +
				geom_abline(intercept=0, slope=1, linetype="solid") +  
				xlim(c(0, max.xy)) +  
				ylim(c(0, max.xy)) +
				xlab("Total predicted events by decile of predicted risk") +
				ylab("Total observed events by decile of predicted risk") + 
				theme(panel.grid.major = element_blank(), 
						axis.line = element_line(colour = "black"), 
						plot.margin = margin(2,2,2,2, "cm"),
						panel.grid.minor = element_blank(), 
						panel.background = element_blank(), 
						text=element_text(size=22)) +
				ggtitle(title.in) 
		return(p)
	}
	
	plot.calibration.hoslem <- function(observed, expected) {
		obsexp <- data.table(observed, expected)
		obsexp[, riskdecile.group := cut(expected,
						quantile(x=expected, probs=seq(0, 1, by=0.1)))]
		expected.riskdecile <- obsexp[, sum(expected), by=riskdecile.group][["V1"]]
		observed.riskdecile <- obsexp[, sum(observed), by=riskdecile.group][["V1"]]
		teststat.chisq <- sum((observed.riskdecile - 
							expected.riskdecile) ^2 / expected.riskdecile)
		df <- 10
		pvalue <- 1 - pchisq(q=teststat.chisq, df=df)
		cat("Hosmer-Lemeshow chi-square statistic", 
				teststat.chisq, "with", df, "degrees of freedom",
				"p = ", signif(pvalue, 1), "\n")
		
		calib <- data.table(x=expected.riskdecile, y=observed.riskdecile)
		max.xy <- max(c(calib$x, calib$y))
		
		p <- ggplot(data=calib, aes(x=x, y=y)) + geom_point() +
				geom_abline(intercept=0, slope=1, linetype="dotted") +  
				xlim(c(0, max.xy)) +  
				ylim(c(0, max.xy)) +
				xlab(stringr::str_wrap(
								"Expected events by decile of predicted risk", 35)
				) +
				ylab("Observed events")
		return(p)
	}
	
	
	
	plot.calibration.cumulative <- function(observed, expected) {
		obsexp <- data.frame(observed, expected)
		obsexp <- obsexp[order(expected), ]
		obsexp$cumulative.exp <- cumsum(obsexp$exp)
		cumulativeexp.group <- cut(obsexp$cumulative.exp,
				sum(obsexp$expected) * seq(0, 1, by=0.1))
		observed.group <- tapply(obsexp$observed, cumulativeexp.group, sum)
		calib <- data.frame(x=1:10, y=observed.group)
		expected.bin <- 0.1 * sum(obsexp$expected)
		p <- ggplot(data=calib, aes(x=x, y=y)) + geom_point() +
				geom_abline(intercept=expected.bin, slope=0, linetype="dotted") +
				geom_ribbon(aes(ymin=qpois(p=0.25, lambda=expected.bin),
								ymax=qpois(p=0.75, lambda=expected.bin)), fill="grey", alpha=0.3) +
				scale_x_continuous(breaks=1:10) +
				xlab(stringr::str_wrap(
								paste("Ranking of predicted risk grouped into 10 bins of",
										round(expected.bin, 1), "expected events"),
								40)) +
				ylab("Observed events")
		return(p)
	}
	
	table.calibration.cumulative <- function(observed, expected) {
		obsexp <- data.frame(observed, expected)
		obsexp <- obsexp[order(expected), ]
		obsexp$cumulative.exp <- cumsum(obsexp$exp)
		cumulativeexp.group <- cut(obsexp$cumulative.exp,
				sum(obsexp$expected) * seq(0, 1, by=0.1))
		observed.group <- tapply(obsexp$observed, cumulativeexp.group, sum)
		calib <- data.frame(x=1:10, y=observed.group)
		expected.bin <- 0.1 * sum(obsexp$expected)
		table.bins <- as.data.frame(rbind(observed.group, rep(expected.bin, 10)))
		colnames(table.bins) <- 1:10
		rownames(table.bins) <- c("Observed", "Expected")
		return(table.bins)
	}
	
	
	
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
		full.cols <- gsub("albuminuric.grade.max", "albuminuric.grade.mm", full.cols)
		full.cols <- gsub("retgrade", "retgrade.mm", full.cols)
	}
	
	# Save models info for documentation
	saveRDS(full.cols, dataPath("full.cols.forward.RData"))
	saveRDS(base.cols, dataPath("base.cols.forward.RData"))
	saveRDS(with(xx <- data.table(full.cols), 
					print(varsToStrings(xx,"full.cols"))), 
			        dataPath("full.cols.forward.txt.RData"))
	
	# Finalise model dataset
	dataset <- na.omit(surv.obj$retcols.regex(
					paste0(paste(full.cols,collapse="|"), 
					"|slice|entry.age|^age.|sex|duration|slice|study|year|", outcome)))
	
	# Retrieve imputed average model matrix data for factors
	if(opts$options()$impute.baseline) {
		simd.mm <- as.model.matrix(na.omit(surv.obj$retcols.regex("simd")), 
				"simd", c('SIMD2','SIMD3','SIMD4','SIMD5'))[,-1]
		albuminuric.grade.mm <- 
				as.model.matrix(na.omit(surv.obj$retcols.regex("albuminuric.grade.max")), 
						"albuminuric.grade.max", c('Micro','Macro'))[,-1]
		retgrade.mm <- as.model.matrix(na.omit(surv.obj$retcols.regex("retgrade")), 
				"retgrade", c('Non Referable', 'Referable / Eye Clinic'))[,-1]
	}
	
	# Narrow the model columns for adjustment - only continuous values
	adj.cols <- unique(c(full.cols[grepl("med|avg|crt", full.cols)],
					base.cols[grepl("med|avg|crt", base.cols)]))
	
	
	# Add squared terms to continuous data and squared logs for skewed data
	var.skewed = ""
	var.remove = as.character()	
	if(opts$options()$check.skewness) {
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
			data = dataset
	)
	saveRDS(base.model, dataPath(paste0("base.model.Rdata")))
	
	
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
			data = dataset
	)
	
	# Required for model environment dataset
	.GlobalEnv$dataset <- dataset
	
	# Run selection
	set.seed(opts$options()$seed)
	forward.model <-
			createOrRetrieve(
					"forward.model.Rdata",
					creator = stepAIC,
					creator.args = list(
							object = base.model,
							scope = list(lower=base.model, upper=full.model),
							direction = "forward",
							trace = 1,
							k = 3,
							method = "approximate"
					), skip.global = TRUE 
			)
	
	# A quick sanity check - we test model C-statistic and calibration
	# against the fitted model. We test the original model fitted
	# to 365 day intervals with offsets and reverse the offsets by multiplying
	# by the proportion of t complete. 
	
	
	vars.regex <- paste0(opts$options()$model.candidates, "|",opts$options()$model.outcome,
			"|", opts$options()$model.min, "|^slice|entry|exit|serialno|date.of.birth")
	
	data.entry <- surv.obj$retcols.regex(vars.regex)[slice==1]
	data.exit <- surv.obj$retcols.regex(vars.regex)[slice.rev==1]
	
	data.surv <- copy(data.entry)
	data.surv[, study.start := slice.start.date]
	data.surv$study.end <-  data.exit$slice.end.date
	data.surv[, c("slice.start.date", "slice.end.date", "slice", "slice.rev") := NULL]
	
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
					"na.action","effects","y","offset","predictors", "qr",
					"linear.predictors","model")] <- NULL
	save(forward.model.redact, file=dataPath("scottish.forward.model.RData"))
	
	# Prediction performance
	dataset <- surv.obj$retcols.regex(
			paste0(paste(full.cols,collapse="|"), 
					"|serialno|age.crt|entry|year|sex|*.med|*.avg|duration|",
					"retgrade|simd|albuminuric|smoker.evr.max|slice|study|",
					"date.of.birth|", outcome)) 
	stopifnot(nrow(dataset) == nrow(surv.obj$surv.table))
	
	
	set.seed(opts$options()$seed)
	fold.indices <- caret::createFolds(dataset[[outcome]], k = folds.n)
	fold.results.test <- fold.results.train <- full.densities <- list()
	fold.auc <- vector()
	
	#' Overload the step program to use the createOrRetrieve
	#' 
	#' 
	stepw  <- function(
			traindata = traindata,
			fold.indices = fold.indices,
			...
	) 
	{ 
		return(step(...));
	}
	
	#' Perform cross fold regression - parallel capable function
	#' with backup loading via createOrRetrieve
	#'
	stepFold <- function(i, traindata, fold.indices) {
		
		fold <- names(fold.indices[i])
		
		# Retrieve imputed average model matrix data for factors
		simd.mm <- as.model.matrix(traindata[,.(simd)], 
				"simd", c('SIMD2','SIMD3','SIMD4','SIMD5'))[,-1]
		albuminuric.grade.mm <- as.model.matrix(traindata[,.(albuminuric.grade.max)], 
				"albuminuric.grade.max", c('Micro','Macro'))[,-1]
		retgrade.mm <- as.model.matrix(traindata[,.(retgrade)], 
				"retgrade", c('Non Referable', 'Referable / Eye Clinic'))[,-1]
		
		# Base model with offset term
		base.model <- glm(
				formula(paste0(outcome, " ~  offset(log(slice.fu)) + ",  
								paste(base.cols, collapse="+"))), 
				family="poisson",
				data = traindata
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
	registerDoMC(cores=3)
	
	fold.model.train <- foreach(i = 1:folds.n) %dopar% {
		
		stepFold(
				i, 
				traindata = dataset[!fold.indices[[i]]], 
				fold.indices
		)
		
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
		
		dataset[, smoker.mm := smoker.evr.max ]
		
		# The predict function automatically calculates an adjusted probability
		# for intervals of length < 1.
		hazard <- predict.glm( model.in, 
				newdata=dataset[fold.indices[[paste0("Fold",nf)]]], type="response")		
		
		# Hazard for full interval
		hazard <- hazard / dataset[fold.indices[[paste0("Fold",nf)]]]$slice.fu
		
		# Training fold prior probability of an event
		prior.prob.calc <- mean(dataset[!fold.indices[[paste0("Fold",nf)]]][[outcome]]) 
		
		# Store the fold used
		dataset[fold.indices[[paste0("Fold",nf)]], fold := nf]
		
		# Store the unadjusted hazard
		dataset[fold.indices[[paste0("Fold",nf)]], hazard.rate.1yr := hazard ]  
		
	}
	
	
	# Calculate the increment in prediction from the base
	training.base.model.folds <- list()
	for(i in 1:folds.n) {
		training.base.model.folds[[i]] <- glm(
				formula(paste0(outcome, " ~ offset(log(slice.fu)) + ",  
								paste(base.cols, collapse="+"))), 
				family="poisson",
				data = na.omit(dataset[-fold.indices[[i]]])
		)
	}
	
	
	# Predictive performance - testing folds 
	foldModbase <- wevid.foldbase <- list()
	for (nf in 1:folds.n) { 
		
		# Calculate predictions across all folds for all testees.
		hazard.base <- predict(training.base.model.folds[[nf]], 
				newdata=dataset[fold.indices[[nf]]],  type="response")
		# Hazard for full interval
		hazard.base.adj <- hazard.base / dataset[fold.indices[[nf]]]$slice.fu
		dataset[fold.indices[[nf]], hazard.base := hazard.base.adj]
		
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
	

	
	rDiabPrintf("Expanding the survival object to reduce offset effects\n")
	dataset[, study.start := slice.start.date]
	dataset[, study.end := slice.end.date]
	dataset[, c("slice.start.date", "slice.end.date", "slice", "slice.rev") := NULL]
	dataset[, serialno.orig := serialno]
	dataset[, serialno := .I]
	
	calib.obj <- Survivor$new(
			surv.table = dataset, study.start = "study.start",
			study.end = "study.end", id = "serialno"
	)
	calib.obj$generate(interval.short)
	
	# Stop outcome being duplicated
	calib.obj$surv.table[, cvd.any.bin := ifelse(slice.rev==1 & cvd.any.bin==1, 1, 0)]
	sum(calib.obj$surv.table$cvd.any.bin)
	
	
	calib.obj$surv.table[, interval.length := 
					as.integer(slice.end.date - slice.start.date + 1) / 365.25]
	calib.obj$surv.table[, expected := hazard.rate.1yr * interval.length]
	calib.obj$surv.table[, expected.base := hazard.base * interval.length]
	

	interval.uncens.length <- interval.short/365.25
	cat("expected with", interval.short, "day intervals ignoring censoring",
			sum(calib.obj$surv.table$hazard.rate.1yr * interval.uncens.length), "\n") 
	
	# Calculate risk over short interval from hazard rate
	# Probability of event before end of interval
	calib.obj$surv.table[, pred.prob := 
					1 - exp(-hazard.rate.1yr * interval.uncens.length)] 
	cat("sum of probabilities with", interval.short, 
			"day intervals ignoring censoring", sum(calib.obj$surv.table$p), "\n")
	calib.obj$surv.table[, pred.base := 
					1 - exp(-hazard.base * interval.uncens.length)] 
	calib.obj$surv.table[, prior.prob := expected]
	
	calib.obj$surv.table[, roc.age.band := cut(age.crt,breaks=c(20,40,60,+Inf))]
	levels(calib.obj$surv.table$roc.age.band) <- c("20-40","40-60","60+")
	surv.obj$surv.table[, roc.age.band := cut(age.crt,breaks=c(20,40,60,+Inf))]
	levels(surv.obj$surv.table$roc.age.band) <- c("20-40","40-60","60+")
	
	
	model.roc.ci <- list()
	model.roc.ci[["all"]] <- pROC::ci.auc(calib.obj$surv.table[[outcome]],
			calib.obj$surv.table$pred.prob)
	model.roc.ci[["male"]] <- pROC::ci.auc(calib.obj$surv.table[sex=='Male'][[outcome]],
			calib.obj$surv.table[sex=='Male']$pred.prob)
	model.roc.ci[["female"]] <- pROC::ci.auc(calib.obj$surv.table[sex=='Female'][[outcome]],
			calib.obj$surv.table[sex=='Female']$pred.prob)
	model.roc.ci[["20-40"]] <- pROC::ci.auc(calib.obj$surv.table[roc.age.band=='20-40'][[outcome]],
			calib.obj$surv.table[roc.age.band=='20-40']$pred.prob)
	model.roc.ci[["40-60"]] <- pROC::ci.auc(calib.obj$surv.table[roc.age.band=='40-60'][[outcome]], 
			calib.obj$surv.table[roc.age.band=='40-60']$pred.prob)
	model.roc.ci[["60+"]] <- pROC::ci.auc(calib.obj$surv.table[roc.age.band=='60+'][[outcome]], 
			calib.obj$surv.table[roc.age.band=='60+']$pred.prob)	
	saveRDS(model.roc.ci, 
			dataPath("confint.forward.auroc.full.Rdata"))
	
	model.roc.ci.base <- list()
	model.roc.ci.base[["all"]] <- pROC::ci.auc(calib.obj$surv.table[[outcome]],
			calib.obj$surv.table$pred.base)
	model.roc.ci.base[["male"]] <- pROC::ci.auc(calib.obj$surv.table[sex=='Male'][[outcome]],
			calib.obj$surv.table[sex=='Male']$pred.base)
	model.roc.ci.base[["female"]] <- pROC::ci.auc(calib.obj$surv.table[sex=='Female'][[outcome]], 
			calib.obj$surv.table[sex=='Female']$pred.base)
	model.roc.ci.base[["20-40"]] <- pROC::ci.auc(calib.obj$surv.table[roc.age.band=='20-40'][[outcome]], 
			calib.obj$surv.table[roc.age.band=='20-40']$pred.base)
	model.roc.ci.base[["40-60"]] <- pROC::ci.auc(calib.obj$surv.table[roc.age.band=='40-60'][[outcome]],
			calib.obj$surv.table[roc.age.band=='40-60']$pred.base)
	model.roc.ci.base[["60+"]] <- pROC::ci.auc(calib.obj$surv.table[roc.age.band=='60+'][[outcome]],
			calib.obj$surv.table[roc.age.band=='60+']$pred.base)	
	saveRDS(model.roc.ci.base, 
			dataPath("confint.forward.auroc.base.Rdata"))
	
	
	wevid.result.1yr.female <- with(calib.obj$surv.table[sex=="Female"],
			Wdensities(
					cvd.any.bin,
					pred.prob,
					prior.prob
			)
	)
	
	wevid.result.1yr.male <- with(calib.obj$surv.table[sex=="Male"],
			Wdensities(
					cvd.any.bin,
					pred.prob,
					prior.prob
			)
	)
	
	
	wevid.result.1yr.20to40 <- with(calib.obj$surv.table[roc.age.band=='20-40'],
			Wdensities(
					cvd.any.bin,
					pred.prob,
					prior.prob
			)
	)
	
	wevid.result.1yr.40to60 <- with(calib.obj$surv.table[roc.age.band=='40-60'],
			Wdensities(
					cvd.any.bin,
					pred.prob,
					prior.prob
			)
	)
	
	wevid.result.1yr.60P <- with(calib.obj$surv.table[roc.age.band=='60+'],
			Wdensities(
					cvd.any.bin,
					pred.prob,
					prior.prob
			)
	)
	
	wevid.result.1yrb.female <- with(calib.obj$surv.table[sex=="Female"],
			Wdensities(
					cvd.any.bin,
					pred.base,
					prior.base
			)
	)
	
	wevid.result.1yrb.male <- with(calib.obj$surv.table[sex=="Male"],
			Wdensities(
					cvd.any.bin,
					pred.base,
					prior.base
			)
	)
	
	wevid.result.1yrb.20to40 <- with(calib.obj$surv.table[roc.age.band=='20-40'],
			Wdensities(
					cvd.any.bin,
					pred.base,
					prior.base
			)
	)
	
	wevid.result.1yrb.40to60 <- with(calib.obj$surv.table[roc.age.band=='40-60'],
			Wdensities(
					cvd.any.bin,
					pred.base,
					prior.base
			)
	)
	
	wevid.result.1yrb.60P <- with(calib.obj$surv.table[roc.age.band=='60+'],
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
					`Crude Î› (bits)` = `Crude Î› (bits).x` - `Crude Î› (bits).y`,
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
			"C-statistic (Final) (95% CI)")
	setnames(all.wevid.side, "Crude C-statistic.y", 
			"C-statistic (Base) (95% CI)")
	setnames(all.wevid.side, "Crude Î› (bits).x", 
			"Î› (Final)")
	setnames(all.wevid.side, "Crude Î› (bits).y", 
			"Î› (Base)")
	
# Format and order columns for paper
	all.wevid.side <- cbind(all.wevid.side, abs(all.wevid.base[,7]-all.wevid.full[,7]))
	fmat <- function(x) {return(format(round(x,2),snmall=2))}
	all.wevid.side <- all.wevid.side[, c(Group=list(Group),lapply(.SD, fmat)), .SDcols=c(2:6)]
	all.wevid.side <- all.wevid.side[,c(1,3,2,5,4,6)]
	setnames(all.wevid.side,"Test log-likelihood (nats)", 
			"log-likelihood increment (natural log units)")
	saveRDS(all.wevid.side, dataPath("all.wevid.side.Rdata"))
	

	cp1 <- (plot.calibration(
			observed=calib.obj$surv.table[[opts$options()$model.outcome]],
			expected=calib.obj$surv.table$pred.prob, 
			title.in = paste0("All (",
			round(sum(surv.obj$surv.table[slice==1]$all.fu)/365.25)," person/yrs)")))
	
	cp2 <- plot.calibration(observed=calib.obj$surv.table[
					sex=='Female'][[opts$options()$model.outcome]],
			expected=calib.obj$surv.table[sex=='Female']$pred.prob, 
			title.in = paste0("Female (",
					round(sum(surv.obj$surv.table[sex=='Female' 
						& slice==1]$all.fu)/365.25) ," person/yrs)"))
	
	cp3 <- plot.calibration(observed=calib.obj$surv.table[
					sex=='Male'][[opts$options()$model.outcome]],
			expected=calib.obj$surv.table[sex=='Male']$pred.prob, 
			title.in = paste0("Male (",
					round(sum(surv.obj$surv.table[sex=='Male' & 
							slice==1]$all.fu)/365.25)," person/yrs)"))
	
	cp4 <- plot.calibration(observed=calib.obj$surv.table[
					roc.age.band=='20-40'][[opts$options()$model.outcome]],
			expected=calib.obj$surv.table[roc.age.band=='20-40']$pred.prob, 
			title.in = paste0("Age at entry 20-40 (",
					round(sum(surv.obj$surv.table[roc.age.band=='20-40' & 
							slice==1]$all.fu)/365.25)," person/yrs)"))
	
	cp5 <- plot.calibration(observed=calib.obj$surv.table[
					roc.age.band=='40-60'][[opts$options()$model.outcome]],
			expected=calib.obj$surv.table[roc.age.band=='40-60']$pred.prob, 
			title.in = paste0("Age at entry 40-60 (",
					round(sum(surv.obj$surv.table[roc.age.band=='40-60' & 
							slice==1]$all.fu)/365.25)," person/yrs)"))
	
	
	cp6 <- plot.calibration(observed=calib.obj$surv.table[
					roc.age.band=='60+'][[opts$options()$model.outcome]],
			expected=calib.obj$surv.table[roc.age.band=='60+']$pred.prob, 
			title.in = paste0("Age at entry 60+ (",
					round(sum(surv.obj$surv.table[roc.age.band=='60+' & 
							slice==1]$all.fu)/365.25)," person/yrs)"))
	
	calplots <- list(
			cp1 = cp1,
			cp2 = cp2,
			cp3 = cp3,
			cp4 = cp4,
			cp5 = cp5,
			cp6 = cp6		
	)
	
	saveRDS(calplots, dataPath("calplots.Rdata"))
	
	pdf(file = repPath("calibration_scotland.pdf"), width = 21, height = 24)
	print(ggpubr::ggarrange(
					
					cp1,cp2,cp3,cp4,cp5,cp6,
					ncol = 2, nrow = 3
			
			))
	dev.off()
	
	
	hoslem.forward <- with(calib.obj$surv.table[,.(plotID=.I,
							Observed=cvd.any.bin,Predicted1=pred.prob/slice.fu)], 
			ResourceSelection::hoslem.test(Observed, Predicted1), g=2)
	saveRDS(hoslem.forward, dataPath("hoslem.forward.RData"))
}
