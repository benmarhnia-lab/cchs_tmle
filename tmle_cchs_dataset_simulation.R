##################
## Below are codes to generate mock dataset for analyses in the manuscript
## "Do we need flexible machine-learning algorithms to study long-term health 
## impact of low-concentrations fine particulate matter?"
## Note each variable within the dataset was simulated randomly thus no  
## meaningful association is expected in the dataset
##################
indir <- "" 

#number of subjects to simulate
set.seed(1234)
n <- 10000
nms <- c("PM3yr", "time", "year", "nonacc1", 
         "DHH_SEX", "age", 
         "abor", "ms3c", "airsh", "urban", 
         "immigrant", "yrsincan", 
         "bmi5c", "smk6c", "alc6c", "fvg3c", "exc3c",
         "vismin",   
         "lbs3c", "edu4c",
         "incqu_new_imp", "csize", "depend", "depriv", "ethcon", "instab")
bar <- numeric()
for (i in 1:n) {
  lnth <- sample(1:11, 1, prob = c(rep(0.02, 10), 0.8))
  new <- data.frame(uniqid = rep(i, lnth), 
                    c = 1, ## no censoring in the dataset
                    setNames(replicate(length(nms), NA, simplify = F), nms))
  if (lnth==11) {
    new$nonacc1 <- c(rep(0, lnth-1), sample(0:1, 1, prob = c(0.98, 0.02)))
  } else {
    new$nonacc1 <- c(rep(0, lnth-1), 1)
  }
  new$time <- 0:(lnth-1)
  new$year <- new$time + 2005
  new$PM3yr <- rnorm(lnth, mean = 6.4, sd = 2.3)
  new$PM3yr <- ifelse(new$PM3yr < 1.8, 1.8, new$PM3yr)
  
  ## time-varying covariates
  new$age <- rep(sample(18:79, 1), lnth)
  new$csize <- sample(c(1:6), lnth, replace=TRUE)
  for (j in c("depend", "depriv", "ethcon", "instab", "incqu_new_imp")) {
    new[, j] <- sample(c(1:5), lnth, replace=TRUE)
  }
  
  ## time-fixed excoriates
  for (j in c("DHH_SEX", "immigrant", "vismin", "abor")) {
    new[, j] <- rep(sample(c(0, 1), 1), lnth)
  }
  new$yrsincan <- rep(sample(0:50, 1), lnth)* new$immigrant
  for (j in c("fvg3c", "exc3c", "ms3c", "lbs3c")) {
    new[, j] <- rep(sample(1:3, 1), lnth)
  }
  for (j in c("edu4c")) {
    new[, j] <- rep(sample(1:4, 1), lnth)
  }
  for (j in c("bmi5c", "urban")) {
    new[, j] <- rep(sample(1:5, 1), lnth)
  }
  for (j in c("smk6c", "alc6c", "airsh")) {
    new[, j] <- rep(sample(1:6, 1), lnth)
  }
  bar <- rbind(bar, new)
}
write.csv(bar, file.path(indir, "cchs.csv"), row.names = FALSE)
