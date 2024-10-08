############
## Set stage
############
## Load packages
library(data.table)
library(dplyr)
library(SuperLearner) ## Super Leaner wrapper
library(lmtp) ## Longitudinal modifiable treatment policy (Diaz et al. 2021)
library(earth) ## multiple additive regression splines used by Nugent et al. 2023 and Rudolph et al. 2022
library(glmnet) ## Elastic net regression, including lasso, ridge, and elastic net regularized regression
library(xgboost) ## Extreme Gradient Boosting package
library(ranger) ## fast implementation of Random Forest (Breiman 2001)

### dataset formatting and stage steup
############
### setup
indir <- "" 
outdir <- ""

wave <- "1.1" 
time_points <- 10 ## 10-year cumulative incidence rate

### read in the data
dt <- fread(file.path(indir, "cchs.csv"))

## continuous exposure
expsr <- "PM3yr" 
## use log(PM)
dt$trt <- log(dt[, expsr, with=FALSE])

## data cleaning specific to lmtp
dt$time <- dt$time + 1 ## time should start with 1

## Only keep relevant rows 
dt <- dt[time < time_points + 1, ]

## transform categorical variables into factors
ft.vb <- c("bmi5c", "smk6c", "alc6c", "fvg3c", "exc3c", "lbs3c", "edu4c",
           "incqu_new_imp", 
           "abor", "ms3c", "airsh", "urban",
           "csize", "depend", "depriv", "ethcon", "instab")
dt[, (ft.vb) := lapply(.SD, as.factor), .SDcol = ft.vb]

## specify time-fixed baseline variables
basecovs <- c("yrsincan", "age",
              "abor", "ms3c", "airsh", "urban",
              "bmi5c", "smk6c", "alc6c", "fvg3c", "exc3c",
              "lbs3c", "edu4c",
              "DHH_SEX", "immigrant", "vismin" # these are binary variables
)

## specify time-varying varaibles
timecovs <- c("incqu_new_imp", "csize", "depend", "depriv", "ethcon", "instab", 
              "trt")

## convert from long to wide format
dt.w <- dcast(dt, reformulate("time", response = paste0(c("uniqid"), collapse = "+")),
              value.var = c(timecovs, "nonacc1", "c"), sep=".")
dt.w <- left_join(dt.w, subset(dt[dt$year==2005,], select=paste0(c("uniqid", basecovs))),  by = c('uniqid' = 'uniqid'))

## carry outcome forward if the outcome occured before the end of follow-up
dt.w <- event_locf(dt.w, paste0("nonacc1.", 1:time_points))

## reorder the variables and change back to data.frame for the lmtp package
dt.n <- dt.w[, c("uniqid", basecovs, paste(c(timecovs, "c", "nonacc1"),
                                           rep(1:time_points, each=length(timecovs) + 2),
                                           sep=".")), with=FALSE]
dt.n <- as.data.frame(dt.n)
############ 

### main analyses
############
interventions <- c("threhold8.8")
out <- bar <- foo <- numeric()

for (j in 1:length(interventions)) {
  
  int.val <- log(as.numeric(gsub("threhold", "", interventions[j]))) # this is used by s0, s1, s2, s3, and s4.2  
  
  threshold_lmtp <- function(data, trt) {
    a <- data[[trt]]
    int.val * (a - int.val > 0) + a * (a - int.val <= 0)
  }
  
  end_i <- time_points 
  
  A <- paste0("trt.", 1:end_i)
  Y <- paste0("nonacc1.", 1:end_i) 
  
  timecovs <- timecovs[1:6]
  L <- lapply(1:end_i, function(a) paste0(timecovs, ".", a))
  C <- paste0("c.", 1:end_i)
  W <- basecovs ## time-fixed baseline variables
  
  cat("variables used for each time point:", "\n")
  print(create_node_list(A, time_points, L, W, k = 1))
  
  
  ########################
  ####### s1::NICE parametric g-computation (gfoRmula) with log(pm)
  ####### same R code as in https://github.com/benmarhnia-lab/cchs_g_computation
  ########################
  ### three points
  ### 1. data should be in long format and data.table (use dt, not dt.n; setDT(dt))
  ### 2. time must start with 0 (create a new time variable = dt$time -1 for dt)
  ### 3. for gfoRmula version 1.0.2 (released 2/27/2022), the ID variable must be named "id" (names(dt)[1] <- "id")
  
  ########################
  ####### s2::lmtp_tmle with log(pm) and SL.glm only
  ########################
  set.seed(1234)
  
  learners_outcome <- "SL.glm" ## SL.glm is a main-term-only linear regression (no interactions)
  learners_trt <- "SL.glm" ## SL.glm is a main-term-only linear regression (no interactions)
  
  tmle.int <- lmtp_tmle(data = dt.n, trt = A, outcome = Y, baseline = W, time_vary = L, ## data and variables
                        cens = C, 
                        folds = 10, ## The number of folds to be used for cross-fitting
                        k = 1, ## number of historical measurements considered
                        shift = threshold_lmtp, ## A two argument function that specifies how treatment variables should be shifted.
                        mtp = TRUE, ## Is the intervention of interest a modified treatment policy? Default is FALSE. If treatment variables are continuous this should be TRUE.
                        outcome_type = "survival",
                        learners_outcome = learners_outcome,
                        learners_trt = learners_trt
  )
  saveRDS(tmle.int, file.path(outdir, paste0("tmle_cchs_cycle_", wave, "_", interventions[j], "_s2", ".rds")))
  
  print(tmle.int)
  foo <- data.frame(tidy(tmle.int))
  foo$int <- interventions[j]
  foo$time_points <- end_i
  out <- rbind(out, foo)
  foo <- numeric()
  tmle.ref <- lmtp_tmle(data = dt.n, trt = A, outcome = Y, baseline = W, time_vary = L, ## data and variables
                        cens = C,
                        folds = 10, ## The number of folds to be used for cross-fitting
                        k = 1, ## number of historical measurements considered
                        shift = NULL, ## A two argument function that specifies how treatment variables should be shifted.
                        # mtp = TRUE, ## Is the intervention of interest a modified treatment policy? Default is FALSE. If treatment variables are continuous this should be TRUE.
                        outcome_type = "survival",
                        learners_outcome = learners_outcome,
                        learners_trt = learners_trt
  )
  saveRDS(tmle.ref, file.path(outdir, paste0("tmle_cchs_cycle_", wave, "_nc_s2", ".rds")))
  
  print(tmle.ref)
  foo <- data.frame(tidy(tmle.ref))
  foo$int <- "ref"
  foo$time_points <- end_i
  out <- rbind(out, foo)
  foo <- numeric()
  
  print(lmtp_contrast(tmle.ref, ref=tmle.int))
  foo <- lmtp_contrast(tmle.ref, ref=tmle.int)$vals
  foo$int <- interventions[j]
  foo$time_points <- end_i
  foo$method <- "lmtp_tmle with log(pm) and SL.glm only"
  bar <- rbind(bar, foo)
  foo <- numeric()
  
  ########################
  ####### s3::lmtp_tmle with log(PM) and six algorithms
  ########################
  
  learners_outcome <- c("SL.glm", ## SL.glm is a main-term-only linear regression (no interactions)
                        "SL.xgboost", ## Extreme Gradient Boosting package
                        "SL.ranger", ## fast implementation of Random Forest (Breiman 2001)
                        "SL.earth", ## multiple additive regression splines used by Nugent et al. 2023 and Rudolph et al. 2022
                        "SL.mean", ## simple mean
                        "SL.glmnet" ## Elastic net regression, including lasso, ridge, and elastic net regularized regression
  )
  learners_trt <-  c("SL.glm", ## SL.glm is a main-term-only linear regression (no interactions)
                     "SL.xgboost", ## Extreme Gradient Boosting package
                     "SL.ranger", ## fast implementation of Random Forest (Breiman 2001)
                     "SL.earth", ## multiple additive regression splines used by Nugent et al. 2023 and Rudolph et al. 2022
                     "SL.mean", ## simple mean
                     "SL.glmnet" ## Elastic net regression, including lasso, ridge, and elastic net regularized regression
  )
  
  tmle.int <- lmtp_tmle(data = dt.n, trt = A, outcome = Y, baseline = W, time_vary = L, ## data and variables
                        cens = C, 
                        folds = 10, ## The number of folds to be used for cross-fitting
                        k = 1, ## number of historical measurements considered
                        shift = threshold_lmtp, ## A two argument function that specifies how treatment variables should be shifted.
                        mtp = TRUE, ## Is the intervention of interest a modified treatment policy? Default is FALSE. If treatment variables are continuous this should be TRUE.
                        outcome_type = "survival",
                        learners_outcome = learners_outcome,
                        learners_trt = learners_trt
  )
  saveRDS(tmle.int, file.path(outdir, paste0("tmle_cchs_cycle_", wave, "_", interventions[j], "_s3", ".rds")))
  
  print(tmle.int)
  
  foo <- data.frame(tidy(tmle.int))
  foo$int <- interventions[j]
  foo$time_points <- end_i
  out <- rbind(out, foo)
  foo <- numeric() 
  
  
  time1 <- Sys.time() 
  tmle.ref <- lmtp_tmle(data = dt.n, trt = A, outcome = Y, baseline = W, time_vary = L, ## data and variables
                        cens = C,
                        folds = 10, ## The number of folds to be used for cross-fitting
                        k = 1, ## number of historical measurements considered
                        shift = NULL, ## A two argument function that specifies how treatment variables should be shifted.
                        # mtp = TRUE, ## Is the intervention of interest a modified treatment policy? Default is FALSE. If treatment variables are continuous this should be TRUE.
                        outcome_type = "survival",
                        learners_outcome = learners_outcome,
                        learners_trt = learners_trt
  )
  saveRDS(tmle.ref, file.path(outdir, paste0("tmle_cchs_cycle_", wave, "_nc_s3", ".rds")))
  
  print(tmle.ref)
  
  foo <- data.frame(tidy(tmle.ref))
  foo$int <- "ref"
  foo$time_points <- end_i
  out <- rbind(out, foo)
  foo <- numeric()
  print(lmtp_contrast(tmle.ref, ref=tmle.int))
  
  foo <- lmtp_contrast(tmle.ref, ref=tmle.int)$vals
  foo$int <- interventions[j]
  foo$time_points <- end_i
  foo$method <- "lmtp_tmle with log PM and six algorithms"
  bar <- rbind(bar, foo)
  foo <- numeric()
  
}

############