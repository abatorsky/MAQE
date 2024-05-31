
#Purpose of this script:
#Create simulated data for 20,000-participant test set which emulates the target patient population

library(boot)

set.seed(99)

# set sample size
nn=20000

subject <- rep(1:nn)

#create data frame
test.set <- as.data.frame(cbind(subject))

#create covariates
# x1 will represent age - it is a baseline demographic variable that will be constant during the study
test.set$x11 <- rnorm(nn, mean=agemean, sd=agesd)
test.set$x11 <- ifelse(test.set$x11 < 18, 18, test.set$x11)

scalex11 <- mean(test.set$x11)
dividex11 <- sd(test.set$x11)

# scale by mean and sd of test set x11
test.set$x11 <- (test.set$x11 - scalex11)/dividex11

# x2 will represent binary opioid use variable which can change between baseline and intermediate time point
test.set$x21 <-rbinom(nn, 1, x21osp)

# x3 will represent dichotomized promis depression t score (0=no symptoms, 1=symptoms) which may change between baseline and intermediate time point
test.set$x31 <- rbinom(nn, 1, x31osp)

#set column for z
test.set$z <- rnorm(nn,mean=0,sd=1)

#create stage 1 treatment assignments

# A1 is randomly assigned 50/50 

test.set$a1 <- rbinom(nn, 1, pp)
test.set$pa1 <- pp

# x2 will represent binary opioid use variable which can change between baseline and intermediate time point: x22 = stage 2
test.set$x22_int <- inv.logit(x22b0 + x22b1*test.set$x21 + x22b2*test.set$a1)
test.set$x22 <- rbinom(nn, 1, test.set$x22_int)
test.set$x22_int <- NULL

# x3 will represent dichotomized promis depression t score (0=no symptoms, 1=symptoms) which may change between baseline and intermediate time point: x32 = stage 2
test.set$x32_int <- inv.logit(x32b0 + x32b1*test.set$x31 + x32b2*test.set$a1)
test.set$x32 <- rbinom(nn, 1, test.set$x32_int)
test.set$x32_int <- NULL

#create intermediate outcome which determines responder status
#create error term
e1 <- rnorm(nn, 0, e1ss)

test.set$y1 <- beta_0 + beta_1*test.set$x11 + beta_2*test.set$x21 + beta_3*test.set$x31 + beta_4*test.set$a1 + beta_5*test.set$a1*test.set$x21 + beta_6*test.set$a1*test.set$x31 + eta*test.set$z*test.set$a1 + e1 + beta11*test.set$x11*test.set$a1 + beta12*(test.set$x11^2)*test.set$a1 + beta13*(test.set$x11^3)*test.set$a1

respval <- quantile(test.set$y1, 0.6)   #calculate respval 

#create binary variable whether they responded to tx
test.set$resp <- ifelse(test.set$y1 > respval,1,0)

#create second stage tx variable
test.set$a2 <- rbinom(nn, 1, qq)
test.set$pa2 <- qq

#create second stage outcome variable
#create error term
e2 <- rnorm(nn, 0, e2ss)
test.set$y2 <- gma_0 + gma_1*test.set$x11 + gma_2*test.set$x22 + gma_3*test.set$x32 + gma_4*test.set$a2 + gma_5*test.set$x22*test.set$a2 + gma_6*test.set$x32*test.set$a2 + gma_7*test.set$resp + gma_8*test.set$resp*test.set$a2 + gma_9*test.set$a1 + eta*test.set$z*test.set$a2 + e2 + gma11*test.set$x11*test.set$a2 + gma12*(test.set$x11^2)*test.set$a2 + gma13*(test.set$x11^3)*test.set$a2

#add vector of 1's to X data frame
test.set$intercept <- 1

# add noise variables
test.set$noise1 <- rnorm(nn, mean=n1mean, sd=n1sd)
test.set$noise2 <- rnorm(nn, mean=n2mean, sd=n2sd)
test.set$noise3 <- rnorm(nn, mean=n3mean, sd=n3sd)
test.set$noise4 <- rbinom(nn, 1, n4mean)
test.set$noise5 <- rbinom(nn, 1, n5mean)
test.set$noise6 <- runif(nn, min=n6min, max=n6max) + test.set$x11*0.1
test.set$noise7 <- rnorm(nn, mean=n7mean, sd=n7sd)
test.set$noise8 <- rbinom(nn, 1, n8mean)
test.set$noise9 <- rbinom(nn, 1, n9mean)
test.set$noise10 <- rnorm(nn, mean=n10mean, sd=n10sd) + test.set$noise7*0.5

test.set$noise11 <- rnorm(nn, mean=0, sd=1)
test.set$noise12 <- rnorm(nn, mean=0, sd=1)
test.set$noise13 <- rnorm(nn, mean=0, sd=1)
test.set$noise14 <- rnorm(nn, mean=0, sd=1)
test.set$noise15 <- rnorm(nn, mean=0, sd=1)
test.set$noise16 <- rnorm(nn, mean=0, sd=1)
test.set$noise17 <- rnorm(nn, mean=0, sd=1)
test.set$noise18 <- rnorm(nn, mean=0, sd=1)
test.set$noise19 <- rnorm(nn, mean=0, sd=1)
test.set$noise20 <- rnorm(nn, mean=0, sd=1)


##### calculate known optimal tx sequence based on outcome model parameters


# create counterfactual a1 and a2

test.set$a1c <- ifelse(test.set$a1==1, 0, 1)
test.set$a2c <- ifelse(test.set$a2==1, 0, 1)

# create counterfactual y1 outcome data with counterfactual a1 

test.set$y1c <- beta_0 + beta_1*test.set$x11 + beta_2*test.set$x21 + beta_3*test.set$x31 + beta_4*test.set$a1c + beta_5*test.set$a1c*test.set$x21 + beta_6*test.set$a1c*test.set$x31 + eta*test.set$z*test.set$a1c + e1  + beta11*test.set$x11*test.set$a1c + beta12*(test.set$x11^2)*test.set$a1c + beta13*(test.set$x11^3)*test.set$a1c

# create counterfactual stage 2 variables respc, x22c and x32c

test.set$respc <- ifelse(test.set$y1c > respval,1,0)

test.set$x22_int <- inv.logit(x22b0 + x22b1*test.set$x21 + x22b2*test.set$a1c)
test.set$x22c <- rbinom(nn, 1, test.set$x22_int)
test.set$x22_int <- NULL

test.set$x32_int <- inv.logit(x32b0 + x32b1*test.set$x31 + x32b2*test.set$a1c)
test.set$x32c <- rbinom(nn, 1, test.set$x32_int)
test.set$x32_int <- NULL

# create counterfactual y2 outcome data with counterfactual a2 assuming real a1

test.set$y2rc <- gma_0 + gma_1*test.set$x11 + gma_2*test.set$x22 + gma_3*test.set$x32 + gma_4*test.set$a2c + gma_5*test.set$x22*test.set$a2c + gma_6*test.set$x32*test.set$a2c + gma_7*test.set$resp + gma_8*test.set$resp*test.set$a2c + gma_9*test.set$a1 + eta*test.set$z*test.set$a2c + e2 + gma11*test.set$x11*test.set$a2c + gma12*(test.set$x11^2)*test.set$a2c + gma13*(test.set$x11^3)*test.set$a2c


# create counterfactual y2 outcome data with counterfactual a2 assuming counterfactual a1

test.set$y2cc <- gma_0 + gma_1*test.set$x11 + gma_2*test.set$x22c + gma_3*test.set$x32c + gma_4*test.set$a2c + gma_5*test.set$x22c*test.set$a2c + gma_6*test.set$x32c*test.set$a2c + gma_7*test.set$respc + gma_8*test.set$respc*test.set$a2c + gma_9*test.set$a1c + eta*test.set$z*test.set$a2c  + e2 + gma11*test.set$x11*test.set$a2c + gma12*(test.set$x11^2)*test.set$a2c + gma13*(test.set$x11^3)*test.set$a2c


# create y2 outcome data based on counterfactual a1 data but real a2

test.set$y2cr <- gma_0 + gma_1*test.set$x11 + gma_2*test.set$x22c + gma_3*test.set$x32c + gma_4*test.set$a2 + gma_5*test.set$x22c*test.set$a2 + gma_6*test.set$x32c*test.set$a2 + gma_7*test.set$respc + gma_8*test.set$respc*test.set$a2 + gma_9*test.set$a1c + eta*test.set$z*test.set$a2 + e2 + gma11*test.set$x11*test.set$a2 + gma12*(test.set$x11^2)*test.set$a2 + gma13*(test.set$x11^3)*test.set$a2

# create 4 values of y1 + y2 for each scenario

test.set$yrr <- test.set$y1 + test.set$y2
test.set$yrc <- test.set$y1 + test.set$y2rc
test.set$ycc <- test.set$y1c + test.set$y2cc
test.set$ycr <- test.set$y1c + test.set$y2cr

# determine max y value

ydata <- as.data.frame(cbind(test.set$yrr, test.set$yrc, test.set$ycc, test.set$ycr))
names(ydata) <- c("yrr", "yrc", "ycc", "ycr")

test.set$ymax <- colnames(ydata)[apply(ydata,1,which.max)]

# assign known optimal tx sequence based on ymax

test.set$stage1.optimal <- rep(3, nn)
test.set$stage2.optimal <- rep(3, nn)

i=0
for(i in 1:dim(test.set)[1]){
  if(test.set$ymax[i]=="yrr"){
    test.set$stage1.optimal[i] <- test.set$a1[i]
    test.set$stage2.optimal[i] <- test.set$a2[i]
  }else if(test.set$ymax[i]=="yrc"){
    test.set$stage1.optimal[i] <- test.set$a1[i]
    test.set$stage2.optimal[i] <- test.set$a2c[i]
  }else if(test.set$ymax[i]=="ycc"){
    test.set$stage1.optimal[i] <- test.set$a1c[i]
    test.set$stage2.optimal[i] <- test.set$a2c[i]
  }else if(test.set$ymax[i]=="ycr"){
    test.set$stage1.optimal[i] <- test.set$a1c[i]
    test.set$stage2.optimal[i] <- test.set$a2[i]
  }
}

# get maximum value for each subject

test.set <- test.set %>%
  mutate(max_value = pmax(yrr, yrc, ycc, ycr))

# delete unnecessary variables

rm(ydata)

