library(ggplot2)
library(tidyr)
library(boot)
library(dplyr)
library(stringr)

### Purpose of this script:
### Read in simulation parameters from param.csv and use for-loop to run different experimental settings
##### For each set of parameters, source GenerateTestSet.R to create test dataset which represents general population of 20,000 people
##### Determine known optimal treatment sequence based on known outcome model for test set population
##### Calculate value(s) for all 4 one size fits all treatments for test set population
####### Use for-loop to run 500 simulations for each set of parameters
####### Within each loop, create training trial and OS datasets and run different versions of the MAQE and SQE
####### Save summary statistics from 500 simulations including  % correctly classified and value to create box plots

# source required functions and load in parameters
### SET WORKING DIRECTORY HERE ###
setwd("C:/Example")
# load functions
source('Functions.R')

######## Loop: One loop for each experimental parameter set #########

# Load in file with different simulation settings

input.params <- read.csv(file="param.csv")
nExp <- nrow(input.params)

for(param.set in 1:nExp){
  
  # read in parameters for this experiment(loop)
  exp.param <- input.params[param.set,]
  
  exp.name <- exp.param$exp
  x22b0 <- exp.param$x22b0
  x22b1 <- exp.param$x22b1
  x22b2 <- exp.param$x22b2
  x32b0 <- exp.param$x32b0
  x32b1 <- exp.param$x32b1
  x32b2 <- exp.param$x32b2
  pp <- exp.param$pp
  qq <- exp.param$qq
  e1ss <- exp.param$e1ss
  e2ss <- exp.param$e2ss
  beta_0 <- exp.param$beta_0
  beta_1 <- exp.param$beta_1
  beta_2 <- exp.param$beta_2
  beta_3 <- exp.param$beta_3
  beta_4 <- exp.param$beta_4
  beta_5 <- exp.param$beta_5
  beta_6 <- exp.param$beta_6
  beta11 <- exp.param$beta_11
  beta12 <- exp.param$beta_12
  beta13 <- exp.param$beta_13
  beta_m1 <- exp.param$beta_m1
  beta_m2 <- exp.param$beta_m2
  beta_m3 <- exp.param$beta_m3
  gma_0 <- exp.param$gma_0
  gma_1 <- exp.param$gma_1
  gma_2 <- exp.param$gma_2
  gma_3 <- exp.param$gma_3
  gma_4 <- exp.param$gma_4
  gma_5 <- exp.param$gma_5
  gma_6 <- exp.param$gma_6
  gma_7 <- exp.param$gma_7
  gma_8 <- exp.param$gma_8
  gma_9 <- exp.param$gma_9
  gma11 <- exp.param$gma_11
  gma12 <- exp.param$gma_12
  gma13 <- exp.param$gma_13
  gma_m1 <- exp.param$gma_m1
  gma_m2 <- exp.param$gma_m2
  gma_m3 <- exp.param$gma_m3
  eta <- exp.param$eta
  n <- exp.param$n
  m <- exp.param$m
  agemean <- exp.param$agemean
  agesd <- exp.param$agesd
  psb0 <- exp.param$psb0
  psb1 <- exp.param$psb1
  psb2 <- exp.param$psb2
  psb3 <- exp.param$psb3
  rctagemean <- exp.param$rctagemean
  rctagesd <- exp.param$rctagesd
  x21rctp <- exp.param$x21rctp
  x31rctp <- exp.param$x31rctp
  x21osp <- exp.param$x21osp
  x31osp <- exp.param$x31osp
  n1mean <- exp.param$n1mean
  n1sd <- exp.param$n1sd
  n2mean <- exp.param$n2mean
  n2sd <- exp.param$n2sd
  n3mean <- exp.param$n3mean
  n3sd <- exp.param$n3sd
  n4mean <- exp.param$n4mean
  n5mean <- exp.param$n5mean
  n7mean <- exp.param$n7mean
  n7sd <- exp.param$n7sd
  n8mean <- exp.param$n8mean
  n9mean <- exp.param$n9mean
  n6min <- exp.param$n6min
  n6max <- exp.param$n6max
  n10mean <- exp.param$n10mean
  n10sd <- exp.param$n10sd
  omegaw <- exp.param$w
  
  # create output file to save strings of text and results
  textfile.name=paste("Results",formatC(exp.name),".txt",sep="")
  cat(sprintf("Results for experiment %s", exp.name), file=textfile.name, append=F, sep="\n")
  
#########

# source file to create test set of 20,000 people with known optimal treatment sequences based on known outcome model parameters
# random seed is set in GenerateTestSet.R so that the test set is created the same way each time
source('GenerateTestSet.R')

# Calculate values for all 4 one size fits all treatments

# assign indicator for which of the 4 treatment possibilities someone received to compare one size fits all 
# indicator onesize = 1 if tx={0,0}; onesize = 2 if tx={0,1}; onesize = 3 if tx={1,0}; onesize = 4 if tx={1,1}

i=0
for(i in 1:dim(test.set)[1]){
  if(test.set$a1[i]==0 && test.set$a2[i]==0){
    test.set$onesize[i] = 1
  }else if(test.set$a1[i]==0 && test.set$a2[i]==1){
    test.set$onesize[i]=2
  }else if(test.set$a1[i]==1 && test.set$a2[i]==0){
    test.set$onesize[i]=3
  }else if(test.set$a1[i]==1 && test.set$a2[i]==1){
    test.set$onesize[i]=4
  }
}

# calculate value of value function for everyone - remember that pp is prob tx assignment A=1 at stage 1 and qq is prob tx assignment A=1 at stage 2

test.set$value <- (test.set$y1 + test.set$y2)/(test.set$pa1*test.set$pa2)

test.set$vald <- test.set$max_value
oraclevalue <- sum(test.set$vald)/nn

# number of simulations

simnum <- 500

# create empty lists
# MAQE.w0 = MAQE2 in the paper; MAQE.wTP = MAQE1 in the paper; SQE.RCT = SQE3 in the paper; SQE.OS = SQE1 in the paper; SQE.RCTOS = SQE2 in the paper
# note that this script is set up to evaluate the "benefit," an evaluation metric that is not well-defined in the literature but sometimes assumed to be the difference in value 
  #   if everyone followed the estimated optimal treatment sequence vs if everyone followed the opposite. This is a good metric to evaluate ITRs (single-stage treatment rules) but may not
  #   be the best evaluation metric for multi-stage trials because there are many people for whom the estimated optimal DTR only gets part of the optimal sequence right. 

AugMethodVar <- vector("list", 5)
names(AugMethodVar) <- c("MAQE.w0", "MAQE.wTP", "SQE.RCT", "SQE.OS", "SQE.RCTOS")

values <- vector("list", 5)
names(values) <- c("MAQE.w0", "MAQE.wTP", "SQE.RCT", "SQE.OS", "SQE.RCTOS")
values$MAQE.w0 <- data.frame(matrix(vector(), simnum, 2, dimnames=list(c(), c("C1","C2"))),stringsAsFactors=F)
names(values$MAQE.w0) <- c("Value", "Benefit")
values$MAQE.wTP <- data.frame(matrix(vector(), simnum, 2, dimnames=list(c(), c("C1","C2"))),stringsAsFactors=F)
names(values$MAQE.wTP) <- c("Value", "Benefit")
values$SQE.RCT <- data.frame(matrix(vector(), simnum, 2, dimnames=list(c(), c("C1","C2"))),stringsAsFactors=F)
names(values$SQE.RCT) <- c("Value", "Benefit")
values$SQE.OS <- data.frame(matrix(vector(), simnum, 2, dimnames=list(c(), c("C1","C2"))),stringsAsFactors=F)
names(values$SQE.OS) <- c("Value", "Benefit")
values$SQE.RCTOS <- data.frame(matrix(vector(), simnum, 2, dimnames=list(c(), c("C1","C2"))),stringsAsFactors=F)
names(values$SQE.RCTOS) <- c("Value", "Benefit")


seed = 86

for(ii in 1:simnum){
  # create training dataset
  set.seed(seed)
  
  #Simulate SMART
  
  subject <- rep(1:n)
  
  #create data frame
  rct <- as.data.frame(cbind(subject))
  
  #create covariates
  # x1 will represent age - it is a baseline demographic variable that will be constant during the study
  rct$x11 <- rnorm(n, mean=rctagemean, sd=rctagesd)
  # all participants need to be 18+
  rct$x11 <- ifelse(rct$x11 < 18, 18, rct$x11)
  
  #scale by mean of x11 in test set, divide by sd of x11 in test set
  rct$x11 <- (rct$x11 - scalex11)/dividex11
  
  # x2 will represent binary opioid use variable which can change between baseline and intermediate time point
  rct$x21 <- rbinom(n, 1, x21rctp)
  
  # x3 will represent dichotomized promis depression t score (0=no symptoms, 1=symptoms) which may change between baseline and intermediate time point
  rct$x31 <- rbinom(n, 1, x31rctp)
  
  #set column for z
  rct$z <- rnorm(n,mean=0,sd=1)
  
  #create stage 1 treatment assignments
  
  rct$a1 <- rbinom(n, 1, pp)
  rct$pa1 <- pp

  # x2 will represent binary opioid use variable which can change between baseline and intermediate time point: x22 = stage 2
  rct$x22_int <- inv.logit(x22b0 + x22b1*rct$x21 + x22b2*rct$a1)
  rct$x22 <- rbinom(n, 1, rct$x22_int)
  rct$x22_int <- NULL
  
  # x3 will represent dichotomized promis depression t score (0=no symptoms, 1=symptoms) which may change between baseline and intermediate time point: x32 = stage 2
  rct$x32_int <- inv.logit(x32b0 + x32b1*rct$x21 + x32b2*rct$a1)
  rct$x32 <- rbinom(n, 1, rct$x32_int)
  rct$x32_int <- NULL
  
  #create intermediate outcome which determines responder status
  #create error term
  e1 <- rnorm(n, 0, e1ss)
  rct$y1 <- beta_0 + beta_1*rct$x11 + beta_2*rct$x21 + beta_3*rct$x31 + beta_4*rct$a1 + beta_5*rct$a1*rct$x21 + beta_6*rct$a1*rct$x31 + eta*rct$z*rct$a1 + e1 + beta11*rct$x11*rct$a1 + beta12*(rct$x11^2)*rct$a1 + beta13*(rct$x11^3)*rct$a1
  
  
  #create binary variable for whether they responded to tx - resp = 1 if y1 > 60th percentile from test set
  rct$resp <- ifelse(rct$y1 > respval,1,0)
  
  #create second stage tx variable 
  rct$a2 <- rbinom(n, 1, qq)
  rct$pa2 <- qq
  
  #create second stage outcome variable
  #create error term
  e2 <- rnorm(n, 0, e2ss)
  rct$y2 <- gma_0 + gma_1*rct$x11 + gma_2*rct$x22 + gma_3*rct$x32 + gma_4*rct$a2 + gma_5*rct$x22*rct$a2 + gma_6*rct$x32*rct$a2 + gma_7*rct$resp + gma_8*rct$resp*rct$a2 + gma_9*rct$a1 + eta*rct$z*rct$a2 + e2  + gma11*rct$x11*rct$a2 + gma12*(rct$x11^2)*rct$a2 + gma13*(rct$x11^3)*rct$a2
  
  #add vector of 1's to X data frame
  rct$intercept <- 1
  
  #add indicator for being in RCT
  rct$smart <- 1
  
  # add noise variables
  rct$noise1 <- rnorm(n, mean=n1mean, sd=n1sd)
  rct$noise2 <- rnorm(n, mean=n2mean, sd=n2sd)
  rct$noise3 <- rnorm(n, mean=n3mean, sd=n3sd)
  rct$noise4 <- rbinom(n, 1, n4mean)
  rct$noise5 <- rbinom(n, 1, n5mean)
  rct$noise6 <- runif(n, min=n6min, max=n6max) + rct$x11*0.1
  rct$noise7 <- rnorm(n, mean=n7mean, sd=n7sd)
  rct$noise8 <- rbinom(n, 1, n8mean)
  rct$noise9 <- rbinom(n, 1, n9mean)
  rct$noise10 <- rnorm(n, mean=n10mean, sd=n10sd) + rct$noise7*0.5
  
  rct$noise11 <- rnorm(n, mean=0, sd=1)
  rct$noise12 <- rnorm(n, mean=0, sd=1)
  rct$noise13 <- rnorm(n, mean=0, sd=1)
  rct$noise14 <- rnorm(n, mean=0, sd=1)
  rct$noise15 <- rnorm(n, mean=0, sd=1)
  rct$noise16 <- rnorm(n, mean=0, sd=1)
  rct$noise17 <- rnorm(n, mean=0, sd=1)
  rct$noise18 <- rnorm(n, mean=0, sd=1)
  rct$noise19 <- rnorm(n, mean=0, sd=1)
  rct$noise20 <- rnorm(n, mean=0, sd=1)

  
  #Simulate Observational Study
  
  ##### Simulate OS #####
  subject <- rep((n+1):(n+m))
  
  #create data frame
  os <- as.data.frame(cbind(subject))
  
  #create covariates
  # x1 will represent age - it is a baseline demographic variable that will be constant during the study
  
  # x1 for OS participants
  os$x11 <- rnorm(m, mean=agemean, sd=agesd)
  # all participants need to be 18+
  os$x11 <- ifelse(os$x11 < 18, 18, os$x11)
  
  #scale by mean of x11 in test set, divide by sd of x11 in test set
  os$x11 <- (os$x11 - scalex11)/dividex11
  
  # x2 will represent binary opioid use variable which can change between baseline and intermediate time point
  os$x21 <-rbinom(m, 1, x21osp)
  
  # x3 will represent dichotomized promis depression t score (0=no symptoms, 1=symptoms) which may change between baseline and intermediate time point
  os$x31 <- rbinom(m, 1, x31osp)
  
  # create variable z, an unmeasured confounder
  os$z <- rnorm(m,mean=0,sd=1)
  
  #create treatment assignments
  # A1 is generated for OS population - use propensity score based on X1-X3
  # use logistic regression parameters created above - apply inverse logit to find probability of y=1, then randomly assign 0 or 1 based on prob
  os$pa1 <- inv.logit(psb0 + psb1*os$x11 + psb2*os$x21 +psb3*os$x31 + 2*os$z)
  os$a1 <- rbinom(m, 1, os$pa1)
  
  # x2 will represent binary opioid use variable which can change between baseline and intermediate time point
  os$x22_int <- inv.logit(x22b0 + x22b1*os$x21 + x22b2*os$a1)
  os$x22 <- rbinom(m, 1, os$x22_int)
  os$x22_int <- NULL
  
  # x3 will represent dichotomized promis depression t score (0=no symptoms, 1=symptoms) which may change between baseline and intermediate time point
  os$x32_int <- inv.logit(x32b0 + x32b1*os$x21 + x32b2*os$a1)
  os$x32 <- rbinom(m, 1, os$x32_int)
  os$x32_int <- NULL
  
  #create intermediate outcome which determines responder status
  e1os <- rnorm(m, 0, e1ss)
  os$y1 <- beta_0 + beta_1*os$x11 + beta_2*os$x21 + beta_3*os$x31 + beta_4*os$a1 + beta_5*os$a1*os$x21 + beta_6*os$a1*os$x31 + eta*os$z*os$a1 + e1os + beta11*os$x11*os$a1 + beta12*(os$x11^2)*os$a1 + beta13*(os$x11^3)*os$a1
  
  #create binary variable for whether they responded to tx - resp = 1 if y1 > 60th percentile from test set
  os$resp <- ifelse(os$y1 > respval,1,0)
  
  #create second stage tx variable
  # A1 is generated for OS population - use propensity score based on X1-X3
  # use logistic regression parameters created above - apply inverse logit to find probability of y=1, then randomly assign 0 or 1 based on prob
  os$pa2 <- inv.logit(psb0 + psb1*os$x11 + psb2*os$x22 +psb3*os$x32 + 2*os$z)
  os$a2 <- rbinom(m, 1, os$pa2)
  
  #create second stage outcome variable
  e2os <- rnorm(m, 0, e2ss)
  os$y2 <- gma_0 + gma_1*os$x11 + gma_2*os$x22 + gma_3*os$x32 + gma_4*os$a2 + gma_5*os$x22*os$a2 + gma_6*os$x32*os$a2 + gma_7*os$resp + gma_8*os$resp*os$a2 + gma_9*os$a1 + eta*os$z*os$a2 + e2os  + gma11*os$x11*os$a2 + gma12*(os$x11^2)*os$a2 + gma13*(os$x11^3)*os$a2
  
  #add vector of 1's to X data frame
  os$intercept <- 1
  
  #add indicator for being in RCT
  os$smart <- 0
  
  # add noise variables
  os$noise1 <- rnorm(m, mean=n1mean, sd=n1sd)
  os$noise2 <- rnorm(m, mean=n2mean, sd=n2sd)
  os$noise3 <- rnorm(m, mean=n3mean, sd=n3sd)
  os$noise4 <- rbinom(m, 1, n4mean)
  os$noise5 <- rbinom(m, 1, n5mean)
  os$noise6 <- runif(m, min=n6min, max=n6max) + os$x11*0.1
  os$noise7 <- rnorm(m, mean=n7mean, sd=n7sd)
  os$noise8 <- rbinom(m, 1, n8mean)
  os$noise9 <- rbinom(m, 1, n9mean)
  os$noise10 <- rnorm(m, mean=n10mean, sd=n10sd) + os$noise7*0.5
  
  os$noise11 <- rnorm(m, mean=0, sd=1)
  os$noise12 <- rnorm(m, mean=0, sd=1)
  os$noise13 <- rnorm(m, mean=0, sd=1)
  os$noise14 <- rnorm(m, mean=0, sd=1)
  os$noise15 <- rnorm(m, mean=0, sd=1)
  os$noise16 <- rnorm(m, mean=0, sd=1)
  os$noise17 <- rnorm(m, mean=0, sd=1)
  os$noise18 <- rnorm(m, mean=0, sd=1)
  os$noise19 <- rnorm(m, mean=0, sd=1)
  os$noise20 <- rnorm(m, mean=0, sd=1)
  
  ########### create one dataset with all data
  
  #######################  USER NOTES ########################## 
  
  # THE ANALYST CAN SPECIFY WHICH DATA TO INCLUDE IN THE MODELS TO ESTIMATE POTENTIAL OUTCOMES
  # THE ANALYST CAN ALSO SPECIFY WHICH DATA TO INCLUDE IN MODEL FOR Q FUNCTION
  
  #######################  USER NOTES ##########################  
  
  simdata <- full_join(rct, os)
  
  # in the param.csv file use a value for the "exp" column that starts with "noise10" to include 10 noise variables in the models for potential outcomes and Q function
  # in the param.csv file use a value for the "exp" column that starts with "noise20" to include 20 noise variables in the models for potential outcomes and Q function
  
  if(startsWith(exp.name, "noise10")){
    mu1var <- list(cbind(simdata$intercept, simdata$x11, simdata$x21, simdata$x31, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10),
                   cbind(simdata$intercept, simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10))
    mu0var <- list(cbind(simdata$intercept, simdata$x11, simdata$x21, simdata$x31, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10),
                   cbind(simdata$intercept, simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10))
    
    Hk <- list(as.data.frame(cbind(simdata$x11, simdata$x21, simdata$x31, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10)),
               as.data.frame(cbind(simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10)))
    MEk <- list(cbind(rct$intercept, rct$x11, rct$x21, rct$x31, rct$noise1, rct$noise2, rct$noise3, rct$noise4, rct$noise5, rct$noise6, rct$noise7, rct$noise8, rct$noise9, rct$noise10),
                cbind(rct$intercept, rct$x11, rct$x22, rct$x32, rct$resp, rct$a1, rct$noise1, rct$noise2, rct$noise3, rct$noise4, rct$noise5, rct$noise6, rct$noise7, rct$noise8, rct$noise9, rct$noise10))

    MEOSk <- list(cbind(os$intercept, os$x11, os$x21, os$x31, os$noise1, os$noise2, os$noise3, os$noise4, os$noise5, os$noise6, os$noise7, os$noise8, os$noise9, os$noise10),
                  cbind(os$intercept, os$x11, os$x22, os$x32, os$resp, os$a1, os$noise1, os$noise2, os$noise3, os$noise4, os$noise5, os$noise6, os$noise7, os$noise8, os$noise9, os$noise10))
    
    simdataMEk <- list(cbind(simdata$intercept, simdata$x11, simdata$x21, simdata$x31, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10),
                   cbind(simdata$intercept, simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10))

    nvars=10
  }else if(startsWith(exp.name, "noise20")){
    mu1var <- list(cbind(simdata$intercept, simdata$x11, simdata$x21, simdata$x31, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10,
                         simdata$noise11, simdata$noise12, simdata$noise13, simdata$noise14, simdata$noise15, simdata$noise16, simdata$noise17, simdata$noise18, simdata$noise19, simdata$noise20),
                   cbind(simdata$intercept, simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10,
                         simdata$noise11, simdata$noise12, simdata$noise13, simdata$noise14, simdata$noise15, simdata$noise16, simdata$noise17, simdata$noise18, simdata$noise19, simdata$noise20))
    mu0var <- list(cbind(simdata$intercept, simdata$x11, simdata$x21, simdata$x31, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10,
                         simdata$noise11, simdata$noise12, simdata$noise13, simdata$noise14, simdata$noise15, simdata$noise16, simdata$noise17, simdata$noise18, simdata$noise19, simdata$noise20),
                   cbind(simdata$intercept, simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10,
                         simdata$noise11, simdata$noise12, simdata$noise13, simdata$noise14, simdata$noise15, simdata$noise16, simdata$noise17, simdata$noise18, simdata$noise19, simdata$noise20))
    
    Hk <- list(as.data.frame(cbind(simdata$x11, simdata$x21, simdata$x31, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10,
                                   simdata$noise11, simdata$noise12, simdata$noise13, simdata$noise14, simdata$noise15, simdata$noise16, simdata$noise17, simdata$noise18, simdata$noise19, simdata$noise20)),
               as.data.frame(cbind(simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10,
                                   simdata$noise11, simdata$noise12, simdata$noise13, simdata$noise14, simdata$noise15, simdata$noise16, simdata$noise17, simdata$noise18, simdata$noise19, simdata$noise20)))
    MEk <- list(cbind(rct$intercept, rct$x11, rct$x21, rct$x31, rct$noise1, rct$noise2, rct$noise3, rct$noise4, rct$noise5, rct$noise6, rct$noise7, rct$noise8, rct$noise9, rct$noise10,
                      rct$noise11, rct$noise12, rct$noise13, rct$noise14, rct$noise15, rct$noise16, rct$noise17, rct$noise18, rct$noise19, rct$noise20),
                cbind(rct$intercept, rct$x11, rct$x22, rct$x32, rct$resp, rct$a1, rct$noise1, rct$noise2, rct$noise3, rct$noise4, rct$noise5, rct$noise6, rct$noise7, rct$noise8, rct$noise9, rct$noise10,
                      rct$noise11, rct$noise12, rct$noise13, rct$noise14, rct$noise15, rct$noise16, rct$noise17, rct$noise18, rct$noise19, rct$noise20))
    
    MEOSk <- list(cbind(os$intercept, os$x11, os$x21, os$x31, os$noise1, os$noise2, os$noise3, os$noise4, os$noise5, os$noise6, os$noise7, os$noise8, os$noise9, os$noise10,
                        os$noise11, os$noise12, os$noise13, os$noise14, os$noise15, os$noise16, os$noise17, os$noise18, os$noise19, os$noise20),
                  cbind(os$intercept, os$x11, os$x22, os$x32, os$resp, os$a1, os$noise1, os$noise2, os$noise3, os$noise4, os$noise5, os$noise6, os$noise7, os$noise8, os$noise9, os$noise10,
                        os$noise11, os$noise12, os$noise13, os$noise14, os$noise15, os$noise16, os$noise17, os$noise18, os$noise19, os$noise20))

    simdataMEk <- list(cbind(simdata$intercept, simdata$x11, simdata$x21, simdata$x31, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10,
                         simdata$noise11, simdata$noise12, simdata$noise13, simdata$noise14, simdata$noise15, simdata$noise16, simdata$noise17, simdata$noise18, simdata$noise19, simdata$noise20),
                   cbind(simdata$intercept, simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1, simdata$noise1, simdata$noise2, simdata$noise3, simdata$noise4, simdata$noise5, simdata$noise6, simdata$noise7, simdata$noise8, simdata$noise9, simdata$noise10,
                         simdata$noise11, simdata$noise12, simdata$noise13, simdata$noise14, simdata$noise15, simdata$noise16, simdata$noise17, simdata$noise18, simdata$noise19, simdata$noise20))
    
    nvars=20
  }else {
    mu1var <- list(cbind(simdata$intercept, simdata$x11, simdata$x21, simdata$x31), cbind(simdata$intercept, simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1))
    mu0var <- list(cbind(simdata$intercept, simdata$x11, simdata$x21, simdata$x31), cbind(simdata$intercept, simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1))
    Hk <- list(as.data.frame(cbind(simdata$x11, simdata$x21, simdata$x31)), as.data.frame(cbind(simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1)))
    MEk <- list(cbind(rct$intercept, rct$x11, rct$x21, rct$x31), cbind(rct$intercept, rct$x11, rct$x22, rct$x32, rct$resp, rct$a1))
    
    MEOSk <- list(cbind(os$intercept, os$x11, os$x21, os$x31), cbind(os$intercept, os$x11, os$x22, os$x32, os$resp, os$a1))
    
    simdataMEk <- list(as.data.frame(cbind(simdata$intercept, simdata$x11, simdata$x21, simdata$x31)), as.data.frame(cbind(simdata$intercept, simdata$x11, simdata$x22, simdata$x32, simdata$resp, simdata$a1)))
    
    nvars=0
  }
  
  
  # inputs for MAQE function
  Yk <- list(simdata$y1, simdata$y2)
  txind <- list(simdata$a1, simdata$a2)
  pi1k <- list(simdata$pa1, simdata$pa2)
  pi0k <- list(1 - simdata$pa1, 1 - simdata$pa2)
  rctind <- simdata$smart
  RCTy <- list(rct$y1, rct$y2)
  RCTa <- list(rct$a1, rct$a2)
  
  # Run SQE with OS:
  OSy <- list(os$y1, os$y2)
  OSa <- list(os$a1, os$a2)
  
  ### run MAQE and SQE functions for each simulation of training trial and OS datasets ###
  # MAQE2 uses w equal to the proportion of trial participants, n/(n+m). This is set in the param.csv file but can be changed. 
  MAQE.wTP <- MAQE(mu1var, mu0var, txind, Hk, Yk, pi1k, pi0k, rctind, MEk, RCTy, RCTa, group="RCT", omegaw=omegaw)
  # MAQE1 uses w equal to 0
  MAQE.w0 <- MAQE(mu1var, mu0var, txind, Hk, Yk, pi1k, pi0k, rctind, MEk, RCTy, RCTa, group="RCT", omegaw=0)
  # SQE3 uses only trial data in the SQE
  SQE.RCT <- SQE(MEk, RCTy, RCTa)
  # SQE1 uses only OS data in the SQE
  SQE.OS <- SQE(MEOSk, OSy, OSa)
  # SQE2 uses both trial and OS data in the SQE
  SQE.RCTOS <- SQE(simdataMEk, Yk, txind)
  
  # evaluate method
  # create history matrices (H1, H2) from the test set to evaluate estimated optimal DTRs from training data
  
  if(startsWith(exp.name, "noise10")){
    h2 <- cbind(test.set$intercept, test.set$x11, test.set$x22, test.set$x32, test.set$resp, test.set$a1, test.set$noise1, test.set$noise2, test.set$noise3, test.set$noise4, test.set$noise5, test.set$noise6, test.set$noise7, test.set$noise8, test.set$noise9, test.set$noise10)
    h1 <- as.matrix(cbind(test.set$intercept, test.set$x11, test.set$x21, test.set$x31, test.set$noise1, test.set$noise2, test.set$noise3, test.set$noise4, test.set$noise5, test.set$noise6, test.set$noise7, test.set$noise8, test.set$noise9, test.set$noise10))
    h2var = 10
    h1var = 10
  }else if(startsWith(exp.name, "noise20")){
    h2 <- cbind(test.set$intercept, test.set$x11, test.set$x22, test.set$x32, test.set$resp, test.set$a1, 
                test.set$noise1, test.set$noise2, test.set$noise3, test.set$noise4, test.set$noise5, test.set$noise6, test.set$noise7, test.set$noise8, test.set$noise9, test.set$noise10,
                test.set$noise11, test.set$noise12, test.set$noise13, test.set$noise14, test.set$noise15, test.set$noise16, test.set$noise17, test.set$noise18, test.set$noise19, test.set$noise20)
    h1 <- as.matrix(cbind(test.set$intercept, test.set$x11, test.set$x21, test.set$x31, 
                          test.set$noise1, test.set$noise2, test.set$noise3, test.set$noise4, test.set$noise5, test.set$noise6, test.set$noise7, test.set$noise8, test.set$noise9, test.set$noise10,
                          test.set$noise11, test.set$noise12, test.set$noise13, test.set$noise14, test.set$noise15, test.set$noise16, test.set$noise17, test.set$noise18, test.set$noise19, test.set$noise20))
    h2var = 20
    h1var = 20
  }else{
    h2 <- cbind(test.set$intercept, test.set$x11, test.set$x22, test.set$x32, test.set$resp, test.set$a1)
    h1 <- as.matrix(cbind(test.set$intercept, test.set$x11, test.set$x21, test.set$x31))
    h2var = 0
    h1var = 0
  }
  

  # 1. MAQE.w0
  
  test.set$MAQE.w0stage2.txrule <- h2%*%as.matrix(t(MAQE.w0$`Stage 2 Coefficients`))
  test.set$MAQE.w0stage2.optimal <- ifelse(test.set$MAQE.w0stage2.txrule > 0, 1, 0)
  
  test.set$MAQE.w0stage1.txrule <- h1%*%as.matrix(t(MAQE.w0$`Stage 1 Coefficients`))
  test.set$MAQE.w0stage1.optimal <- ifelse(test.set$MAQE.w0stage1.txrule > 0, 1, 0)
  
  # assign indicator of whether MAQE.w0 got the "right" sequence
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$MAQE.w0stage1.optimal[i]==test.set$stage1.optimal[i] && test.set$MAQE.w0stage2.optimal[i]==test.set$stage2.optimal[i]){
      test.set$MAQE.w0right[i] = 1
    }else if(test.set$MAQE.w0stage1.optimal[i]!=test.set$stage1.optimal[i] || test.set$MAQE.w0stage2.optimal[i]!=test.set$stage2.optimal[i]){
      test.set$MAQE.w0right[i]=0
    }
  }
  
  AugMethodVar$MAQE.w0[ii] <- length(which(test.set$MAQE.w0right==1))/nn
  
  # assign indicator of whether subject followed estimated optimal tx rule, called opt
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$a1[i]==test.set$MAQE.w0stage1.optimal[i] && test.set$a2[i]==test.set$MAQE.w0stage2.optimal[i]){
      test.set$opt[i] = 1
    }else if(test.set$a1[i]!=test.set$MAQE.w0stage1.optimal[i] || test.set$a2[i]!=test.set$MAQE.w0stage2.optimal[i]){
      test.set$opt[i]=0
    }
  }
  
  # assign indicator of whether subject followed opposite of estimated optimal tx rule, called opp
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$a1[i]!=test.set$MAQE.w0stage1.optimal[i] && test.set$a2[i]!=test.set$MAQE.w0stage2.optimal[i]){
      test.set$opp[i] = 1
    }else if(test.set$a1[i]==test.set$MAQE.w0stage1.optimal[i] || test.set$a2[i]==test.set$MAQE.w0stage2.optimal[i]){
      test.set$opp[i]=0
    }
  }
  
  #calculate value and benefit
  values$MAQE.w0$Value[ii] <- sum(test.set$value[test.set$opt==1])/nn
  values$MAQE.w0$Benefit[ii] <- (sum(test.set$value[test.set$opt==1])-sum(test.set$value[test.set$opp==1]))/nn
  
  # 2. MAQE.wTP
  
  test.set$MAQE.wTPstage2.txrule <- h2%*%as.matrix(t(MAQE.wTP$`Stage 2 Coefficients`))
  test.set$MAQE.wTPstage2.optimal <- ifelse(test.set$MAQE.wTPstage2.txrule > 0, 1, 0)
  
  test.set$MAQE.wTPstage1.txrule <- h1%*%as.matrix(t(MAQE.wTP$`Stage 1 Coefficients`))
  test.set$MAQE.wTPstage1.optimal <- ifelse(test.set$MAQE.wTPstage1.txrule > 0, 1, 0)
  
  # assign indicator of whether MAQE.wTP got the "right" sequence
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$MAQE.wTPstage1.optimal[i]==test.set$stage1.optimal[i] && test.set$MAQE.wTPstage2.optimal[i]==test.set$stage2.optimal[i]){
      test.set$MAQE.wTPright[i] = 1
    }else if(test.set$MAQE.wTPstage1.optimal[i]!=test.set$stage1.optimal[i] || test.set$MAQE.wTPstage2.optimal[i]!=test.set$stage2.optimal[i]){
      test.set$MAQE.wTPright[i]=0
    }
  }
  
  AugMethodVar$MAQE.wTP[ii] <- length(which(test.set$MAQE.wTPright==1))/nn
  
  # assign indicator of whether subject followed estimated optimal tx rule, called opt
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$a1[i]==test.set$MAQE.wTPstage1.optimal[i] && test.set$a2[i]==test.set$MAQE.wTPstage2.optimal[i]){
      test.set$opt[i] = 1
    }else if(test.set$a1[i]!=test.set$MAQE.wTPstage1.optimal[i] || test.set$a2[i]!=test.set$MAQE.wTPstage2.optimal[i]){
      test.set$opt[i]=0
    }
  }
  
  # assign indicator of whether subject followed opposite of estimated optimal tx rule, called opp
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$a1[i]!=test.set$MAQE.wTPstage1.optimal[i] && test.set$a2[i]!=test.set$MAQE.wTPstage2.optimal[i]){
      test.set$opp[i] = 1
    }else if(test.set$a1[i]==test.set$MAQE.wTPstage1.optimal[i] || test.set$a2[i]==test.set$MAQE.wTPstage2.optimal[i]){
      test.set$opp[i]=0
    }
  }
  
  #calculate value and benefit
  values$MAQE.wTP$Value[ii] <- sum(test.set$value[test.set$opt==1])/nn
  values$MAQE.wTP$Benefit[ii] <- (sum(test.set$value[test.set$opt==1])-sum(test.set$value[test.set$opp==1]))/nn
  
  
  # 3. SQE.RCT
  
  test.set$SQE.RCTstage2.txrule <- h2%*%as.matrix(t(SQE.RCT$`Stage 2 Coefficients`))
  test.set$SQE.RCTstage2.optimal <- ifelse(test.set$SQE.RCTstage2.txrule > 0, 1, 0)
  
  test.set$SQE.RCTstage1.txrule <- h1%*%as.matrix(t(SQE.RCT$`Stage 1 Coefficients`))
  test.set$SQE.RCTstage1.optimal <- ifelse(test.set$SQE.RCTstage1.txrule > 0, 1, 0)
  
  # assign indicator of whether SQE.RCT got the "right" sequence
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$SQE.RCTstage1.optimal[i]==test.set$stage1.optimal[i] && test.set$SQE.RCTstage2.optimal[i]==test.set$stage2.optimal[i]){
      test.set$SQE.RCTright[i] = 1
    }else if(test.set$SQE.RCTstage1.optimal[i]!=test.set$stage1.optimal[i] || test.set$SQE.RCTstage2.optimal[i]!=test.set$stage2.optimal[i]){
      test.set$SQE.RCTright[i]=0
    }
  }
  
  AugMethodVar$SQE.RCT[ii] <- length(which(test.set$SQE.RCTright==1))/nn
  
  # assign indicator of whether subject followed estimated optimal tx rule, called opt
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$a1[i]==test.set$SQE.RCTstage1.optimal[i] && test.set$a2[i]==test.set$SQE.RCTstage2.optimal[i]){
      test.set$opt[i] = 1
    }else if(test.set$a1[i]!=test.set$SQE.RCTstage1.optimal[i] || test.set$a2[i]!=test.set$SQE.RCTstage2.optimal[i]){
      test.set$opt[i]=0
    }
  }
  
  # assign indicator of whether subject followed opposite of estimated optimal tx rule, called opp
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$a1[i]!=test.set$SQE.RCTstage1.optimal[i] && test.set$a2[i]!=test.set$SQE.RCTstage2.optimal[i]){
      test.set$opp[i] = 1
    }else if(test.set$a1[i]==test.set$SQE.RCTstage1.optimal[i] || test.set$a2[i]==test.set$SQE.RCTstage2.optimal[i]){
      test.set$opp[i]=0
    }
  }
  
  #calculate value and benefit
  values$SQE.RCT$Value[ii] <- sum(test.set$value[test.set$opt==1])/nn
  values$SQE.RCT$Benefit[ii] <- (sum(test.set$value[test.set$opt==1])-sum(test.set$value[test.set$opp==1]))/nn

  
  # 4. SQE.OS

  test.set$SQE.OSstage2.txrule <- h2%*%as.matrix(t(SQE.OS$`Stage 2 Coefficients`))
  test.set$SQE.OSstage2.optimal <- ifelse(test.set$SQE.OSstage2.txrule > 0, 1, 0)

  test.set$SQE.OSstage1.txrule <- h1%*%as.matrix(t(SQE.OS$`Stage 1 Coefficients`))
  test.set$SQE.OSstage1.optimal <- ifelse(test.set$SQE.OSstage1.txrule > 0, 1, 0)

  # assign indicator of whether SQE.OS got the "right" sequence

  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$SQE.OSstage1.optimal[i]==test.set$stage1.optimal[i] && test.set$SQE.OSstage2.optimal[i]==test.set$stage2.optimal[i]){
      test.set$SQE.OSright[i] = 1
    }else if(test.set$SQE.OSstage1.optimal[i]!=test.set$stage1.optimal[i] || test.set$SQE.OSstage2.optimal[i]!=test.set$stage2.optimal[i]){
      test.set$SQE.OSright[i]=0
    }
  }

  AugMethodVar$SQE.OS[ii] <- length(which(test.set$SQE.OSright==1))/nn

  # assign indicator of whether subject followed estimated optimal tx rule, called opt

  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$a1[i]==test.set$SQE.OSstage1.optimal[i] && test.set$a2[i]==test.set$SQE.OSstage2.optimal[i]){
      test.set$opt[i] = 1
    }else if(test.set$a1[i]!=test.set$SQE.OSstage1.optimal[i] || test.set$a2[i]!=test.set$SQE.OSstage2.optimal[i]){
      test.set$opt[i]=0
    }
  }

  # assign indicator of whether subject followed opposite of estimated optimal tx rule, called opp

  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$a1[i]!=test.set$SQE.OSstage1.optimal[i] && test.set$a2[i]!=test.set$SQE.OSstage2.optimal[i]){
      test.set$opp[i] = 1
    }else if(test.set$a1[i]==test.set$SQE.OSstage1.optimal[i] || test.set$a2[i]==test.set$SQE.OSstage2.optimal[i]){
      test.set$opp[i]=0
    }
  }

  #calculate value and benefit
  values$SQE.OS$Value[ii] <- sum(test.set$value[test.set$opt==1])/nn
  values$SQE.OS$Benefit[ii] <- (sum(test.set$value[test.set$opt==1])-sum(test.set$value[test.set$opp==1]))/nn

  
  # 5. SQE.RCTOS
  
  test.set$SQE.RCTOSstage2.txrule <- h2%*%as.matrix(t(SQE.RCTOS$`Stage 2 Coefficients`))
  test.set$SQE.RCTOSstage2.optimal <- ifelse(test.set$SQE.RCTOSstage2.txrule > 0, 1, 0)
  
  test.set$SQE.RCTOSstage1.txrule <- h1%*%as.matrix(t(SQE.RCTOS$`Stage 1 Coefficients`))
  test.set$SQE.RCTOSstage1.optimal <- ifelse(test.set$SQE.RCTOSstage1.txrule > 0, 1, 0)
  
  # assign indicator of whether SQE.RCTOS got the "right" sequence
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$SQE.RCTOSstage1.optimal[i]==test.set$stage1.optimal[i] && test.set$SQE.RCTOSstage2.optimal[i]==test.set$stage2.optimal[i]){
      test.set$SQE.RCTOSright[i] = 1
    }else if(test.set$SQE.RCTOSstage1.optimal[i]!=test.set$stage1.optimal[i] || test.set$SQE.RCTOSstage2.optimal[i]!=test.set$stage2.optimal[i]){
      test.set$SQE.RCTOSright[i]=0
    }
  }
  
  AugMethodVar$SQE.RCTOS[ii] <- length(which(test.set$SQE.RCTOSright==1))/nn
  
  # assign indicator of whether subject followed estimated optimal tx rule, called opt
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$a1[i]==test.set$SQE.RCTOSstage1.optimal[i] && test.set$a2[i]==test.set$SQE.RCTOSstage2.optimal[i]){
      test.set$opt[i] = 1
    }else if(test.set$a1[i]!=test.set$SQE.RCTOSstage1.optimal[i] || test.set$a2[i]!=test.set$SQE.RCTOSstage2.optimal[i]){
      test.set$opt[i]=0
    }
  }
  
  # assign indicator of whether subject followed opposite of estimated optimal tx rule, called opp
  
  i=0
  for(i in 1:dim(test.set)[1]){
    if(test.set$a1[i]!=test.set$SQE.RCTOSstage1.optimal[i] && test.set$a2[i]!=test.set$SQE.RCTOSstage2.optimal[i]){
      test.set$opp[i] = 1
    }else if(test.set$a1[i]==test.set$SQE.RCTOSstage1.optimal[i] || test.set$a2[i]==test.set$SQE.RCTOSstage2.optimal[i]){
      test.set$opp[i]=0
    }
  }
  
  #calculate value and benefit
  values$SQE.RCTOS$Value[ii] <- sum(test.set$value[test.set$opt==1])/nn
  values$SQE.RCTOS$Benefit[ii] <- (sum(test.set$value[test.set$opt==1])-sum(test.set$value[test.set$opp==1]))/nn
  
  
  seed = seed + 1
}


# save PCC for all methods to create boxplots

boxplotdf <- as.data.frame(cbind(AugMethodVar$MAQE.w0, AugMethodVar$MAQE.wTP, 
                                 AugMethodVar$SQE.RCT, AugMethodVar$SQE.OS, AugMethodVar$SQE.RCTOS))
names(boxplotdf) <- c("MAQE.w0", "MAQE.wTP", "SQE.RCT", "SQE.OS", "SQE.RCTOS")

boxplotdf2 <- boxplotdf %>% gather(Method, PercentCorrect)

# save PCC to create boxplots
CCfile.name=paste("CC",formatC(exp.name),".csv", sep="")
write.csv(boxplotdf, file=CCfile.name, row.names=F)

# save values for all methods to create boxplots

boxplotvalues <- as.data.frame(cbind(values$MAQE.w0$Value, values$MAQE.wTP$Value,
                                     values$SQE.RCT$Value, values$SQE.OS$Value, values$SQE.RCTOS$Value))
names(boxplotvalues) <- c("MAQE.w0", "MAQE.wTP", "SQE.RCT", "SQE.OS", "SQE.RCTOS")

boxplotvalues2 <- boxplotvalues %>% gather(Method, Value)

# save values to create boxplots
valuesfile.name=paste("values",formatC(exp.name),".csv", sep="")
write.csv(boxplotvalues, file=valuesfile.name, row.names=F)

### append summary stats to text file

cat(sprintf("Mean & sd CC of MAQE w=TP %0.3f %0.4f",mean(boxplotdf2$PercentCorrect[boxplotdf2$Method=="MAQE.wTP"]),sd(boxplotdf2$PercentCorrect[boxplotdf2$Method=="MAQE.wTP"])), file=textfile.name, append=T, sep="\n")
cat(sprintf("Mean & sd CC of MAQE w=0 %0.3f %0.4f",mean(boxplotdf2$PercentCorrect[boxplotdf2$Method=="MAQE.w0"]),sd(boxplotdf2$PercentCorrect[boxplotdf2$Method=="MAQE.w0"])), file=textfile.name, append=T, sep="\n")
cat(sprintf("Mean & sd CC of SQE w/RCT %0.3f %0.4f",mean(boxplotdf2$PercentCorrect[boxplotdf2$Method=="SQE.RCT"]),sd(boxplotdf2$PercentCorrect[boxplotdf2$Method=="SQE.RCT"])), file=textfile.name, append=T, sep="\n")
cat(sprintf("Mean & sd CC of SQE.OS %0.3f %0.4f",mean(boxplotdf2$PercentCorrect[boxplotdf2$Method=="SQE.OS"]),sd(boxplotdf2$PercentCorrect[boxplotdf2$Method=="SQE.OS"])), file=textfile.name, append=T, sep="\n")
cat(sprintf("Mean & sd CC of SQE w/RCTOS %0.3f %0.4f",mean(boxplotdf2$PercentCorrect[boxplotdf2$Method=="SQE.RCTOS"]),sd(boxplotdf2$PercentCorrect[boxplotdf2$Method=="SQE.RCTOS"])), file=textfile.name, append=T, sep="\n")

cat(sprintf("Percent of sims MAQE w=TP method performs better (CC) than SQE.RCT: %0.2f",length(which(boxplotdf2$PercentCorrect[boxplotdf2$Method=="MAQE.wTP"] > boxplotdf2$PercentCorrect[boxplotdf2$Method=="SQE.RCT"]))/simnum*100), file=textfile.name, append=T, sep="\n")
cat(sprintf("Percent of sims MAQE w=0 method performs better (CC) than SQE.RCT: %0.2f",length(which(boxplotdf2$PercentCorrect[boxplotdf2$Method=="MAQE.w0"] > boxplotdf2$PercentCorrect[boxplotdf2$Method=="SQE.RCT"]))/simnum*100), file=textfile.name, append=T, sep="\n")
cat(sprintf("Percent of sims MAQE w=TP method performs better (CC) than SQE.RCTOS: %0.2f",length(which(boxplotdf2$PercentCorrect[boxplotdf2$Method=="MAQE.wTP"] > boxplotdf2$PercentCorrect[boxplotdf2$Method=="SQE.RCTOS"]))/simnum*100), file=textfile.name, append=T, sep="\n")
cat(sprintf("Percent of sims MAQE w=0 method performs better (CC) than SQE.RCTOS: %0.2f",length(which(boxplotdf2$PercentCorrect[boxplotdf2$Method=="MAQE.w0"] > boxplotdf2$PercentCorrect[boxplotdf2$Method=="SQE.RCTOS"]))/simnum*100), file=textfile.name, append=T, sep="\n")

cat(sprintf("Value of one size fits all: 0, 0 %0.2f",sum(test.set$value[test.set$onesize==1])/nn), file=textfile.name, append=T, sep="\n")
cat(sprintf("Value of one size fits all: 0, 1 %0.2f",sum(test.set$value[test.set$onesize==2])/nn), file=textfile.name, append=T, sep="\n")
cat(sprintf("Value of one size fits all: 1, 0 %0.2f",sum(test.set$value[test.set$onesize==3])/nn), file=textfile.name, append=T, sep="\n")
cat(sprintf("Value of one size fits all: 1, 1 %0.2f",sum(test.set$value[test.set$onesize==4])/nn), file=textfile.name, append=T, sep="\n")
cat(sprintf("Oracle value: %0.2f", oraclevalue), file=textfile.name, append=T, sep="\n")

cat(sprintf("Mean & sd of Value for MAQE w=TP %0.2f %0.4f",mean(boxplotvalues2$Value[boxplotvalues2$Method=="MAQE.wTP"]),sd(boxplotvalues2$Value[boxplotvalues2$Method=="MAQE.wTP"])), file=textfile.name, append=T, sep="\n")
cat(sprintf("Mean & sd of Value for MAQE w=0 %0.2f %0.4f",mean(boxplotvalues2$Value[boxplotvalues2$Method=="MAQE.w0"]),sd(boxplotvalues2$Value[boxplotvalues2$Method=="MAQE.w0"])), file=textfile.name, append=T, sep="\n")
cat(sprintf("Mean & sd of Value for SQE w/RCT %0.2f %0.4f",mean(boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCT"]),sd(boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCT"])), file=textfile.name, append=T, sep="\n")
cat(sprintf("Mean & sd of Value for SQE.OS %0.2f %0.4f",mean(boxplotvalues2$Value[boxplotvalues2$Method=="SQE.OS"]),sd(boxplotvalues2$Value[boxplotvalues2$Method=="SQE.OS"])), file=textfile.name, append=T, sep="\n")
cat(sprintf("Mean & sd of Value for SQE w/RCTOS %0.2f %0.4f",mean(boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCTOS"]),sd(boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCTOS"])), file=textfile.name, append=T, sep="\n")

cat(sprintf("Percent of sims MAQE w=TP method performs better (value) than best one size fits all tx: %0.2f",length(which(boxplotvalues2$Value[boxplotvalues2$Method=="MAQE.wTP"] > sum(test.set$value[test.set$onesize==1])/nn))/simnum*100), file=textfile.name, append=T, sep="\n")
cat(sprintf("Percent of sims MAQE w=0 method performs better (value) than best one size fits all tx: %0.2f",length(which(boxplotvalues2$Value[boxplotvalues2$Method=="MAQE.w0"] > sum(test.set$value[test.set$onesize==1])/nn))/simnum*100), file=textfile.name, append=T, sep="\n")
cat(sprintf("Percent of sims SQE (w/RCT) method performs better (value) than best one size fits all tx: %0.2f",length(which(boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCT"] > sum(test.set$value[test.set$onesize==1])/nn))/simnum*100), file=textfile.name, append=T, sep="\n")

cat(sprintf("Percent of sims MAQE w=TP method performs better (value) than SQE.RCT: %0.2f",length(which(boxplotvalues2$Value[boxplotvalues2$Method=="MAQE.wTP"] > boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCT"]))/simnum*100), file=textfile.name, append=T, sep="\n")
cat(sprintf("Percent of sims MAQE w=0 method performs better (value) than SQE.RCT: %0.2f",length(which(boxplotvalues2$Value[boxplotvalues2$Method=="MAQE.w0"] > boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCT"]))/simnum*100), file=textfile.name, append=T, sep="\n")
cat(sprintf("Percent of sims MAQE w=TP method performs better (value) than SQE.RCTOS: %0.2f",length(which(boxplotvalues2$Value[boxplotvalues2$Method=="MAQE.wTP"] > boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCTOS"]))/simnum*100), file=textfile.name, append=T, sep="\n")
cat(sprintf("Percent of sims MAQE w=0 method performs better (value) than SQE.RCTOS: %0.2f",length(which(boxplotvalues2$Value[boxplotvalues2$Method=="MAQE.w0"] > boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCTOS"]))/simnum*100), file=textfile.name, append=T, sep="\n")

cat(sprintf("Percent of sims SQE.RCT method performs better (value) than best one size fits all tx: %0.2f",length(which(boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCT"] > sum(test.set$value[test.set$onesize==1])/nn))/simnum*100), file=textfile.name, append=T, sep="\n")
cat(sprintf("Percent of sims SQE.RCTOS method performs better (value) than best one size fits all tx: %0.2f",length(which(boxplotvalues2$Value[boxplotvalues2$Method=="SQE.RCTOS"] > sum(test.set$value[test.set$onesize==1])/nn))/simnum*100), file=textfile.name, append=T, sep="\n")

# add to text file the number of noise variables used and the number of simulations run per setting
cat(c(nvars, h2var, h1var, simnum), file=textfile.name, append=T, sep="\n")

# create file of simulation parameters used for each exp.name as a sanity check
paramlist.name=paste("parameterlist",formatC(exp.name),".csv", sep="")
write.csv(exp.param, file=paramlist.name, row.names=F)

}   # end of loop!
