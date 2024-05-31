# function to perform weighted least squares regression
# takes in covariate matrix (x), covariate matrix weights as a diagonal matrix (w1), outcome weights as a diagonal matrix (w2) and outcome vector (y)
# returns regression coefficients for all variables in x matrix

ls.solver <- function(x, w1, w2, y){
  xt <- t(x)
  xprimewx <- xt%*%w1%*%x
  invxprimewx <- solve(xprimewx)
  xprimey <- xt%*%w2%*%y
  coeff <- invxprimewx%*%xprimey
  return(coeff)
}

# function to create potential outcomes; performs least squares regression, then fits data to the model
# takes in matrices for different treatment groups (mu1var, mu0var), outcome vector (Yk), treatment indicator (txind), trial indicator (trialind), 
#   dataset(s) with which to model potential outcome parameters (group, default="RCT")
# returns a matrix with columns mu1hat and mu0hat, representing the estimated potential outcomes under different treatments

muhats.var <- function(mu1var, mu0var, Yk, txind, trialind, group="RCT"){
  
  dataframe1 <- as.data.frame(cbind(mu1var, Yk, txind, trialind))
  
  if(group=="RCTOS"){
    outcome.model.stage2.a1 <- lm(Yk ~ mu1var[,2:dim(mu1var)[2]], data=dataframe1, subset = txind==1)
    eta.stage2.1 <- as.matrix(outcome.model.stage2.a1$coefficients)
  }else if(group=="OS"){
    outcome.model.stage2.a1 <- lm(Yk ~ mu1var[,2:dim(mu1var)[2]], data=dataframe1, subset = (txind==1 & trialind==0))
    eta.stage2.1 <- as.matrix(outcome.model.stage2.a1$coefficients)
  }else if (group=="RCT"){
    outcome.model.stage2.a1 <- lm(Yk ~ mu1var[,2:dim(mu1var)[2]], data=dataframe1, subset = (txind==1 & trialind==1))
    eta.stage2.1 <- as.matrix(outcome.model.stage2.a1$coefficients)
  }
  
  dataframe0 <- as.data.frame(cbind(mu0var, Yk, txind, trialind))
  
  if(group=="RCTOS"){
    outcome.model.stage2.a0 <- lm(Yk ~ mu0var[,2:dim(mu0var)[2]], data=dataframe0, subset = txind==0)
    eta.stage2.0 <- as.matrix(outcome.model.stage2.a0$coefficients)
  }else if(group=="OS"){
    outcome.model.stage2.a0 <- lm(Yk ~ mu0var[,2:dim(mu0var)[2]], data=dataframe0, subset = (txind==0 & trialind==0))
    eta.stage2.0 <- as.matrix(outcome.model.stage2.a0$coefficients)
  }else if (group=="RCT"){
    outcome.model.stage2.a0 <- lm(Yk ~ mu0var[,2:dim(mu0var)[2]], data=dataframe0, subset = (txind==0 & trialind==1))
    eta.stage2.0 <- as.matrix(outcome.model.stage2.a0$coefficients)
  }
  
  i=0
  mu1.hat <- rep(1:length(Yk))
  mu0.hat <- rep(1:length(Yk))
  for(i in 1:length(Yk)){
    # multiply regression coefficients from model where A=1
    mu1.hat <- as.matrix(mu1var) %*% as.matrix(eta.stage2.1)
    # multiply regression coefficients from model where A=0
    mu0.hat <- as.matrix(mu0var) %*% as.matrix(eta.stage2.0) 
  }
  muhats <-cbind(mu1.hat, mu0.hat) 
  colnames(muhats) <- c("mu1hat","mu0hat")
  return(muhats)
}

# function to run the standard Q-learning estimator
# takes in list of history matrices for multiple stages (SQEHk), list of outcome vectors for multiple stages (SQEYk), list of treatment assignments for multiple stages (SQEAk)
# returns DTR (i.e., list of regression coefficients for each stage corresponding to the covariates from the input)

SQE <- function(SQEHk, SQEYk, SQEAk){
  
  # create weight matrix
  SQEwt2 <- diag(length(SQEYk[[2]]))
  
  # solve for betas stage 2 (tx effects)
  x2ME <- as.matrix(SQEHk[[2]])
  x2TE <- as.matrix(SQEHk[[2]]*SQEAk[[2]])
  x2 <- cbind(x2TE, x2ME)
  SQEcoeff_stage2 <- ls.solver(x2, SQEwt2, SQEwt2, SQEYk[[2]])  # tx effects first coefficients of length of Hk[[2]]
  betas_stage2 <- as.matrix(SQEcoeff_stage2[1:dim(x2ME)[2]]) # stage 2 tx rule
  
  # find optimal stage 2 tx for each subject
  stage2.txrule <- x2ME%*%betas_stage2
  stage2.optimal <- ifelse(stage2.txrule > 0, 1, 0)
  
  # create new dataframe using optimal tx
  SQEHk2df <- as.data.frame(SQEHk[[2]])
  x2optimal <- cbind(SQEHk2df*stage2.optimal, SQEHk2df) # tx effects first, change from 8.7.23
  
  i=0
  stage2Q.hat <- rep(0, length(SQEYk[[2]]))
  for(i in 1:length(SQEYk[[2]])){
    stage2Q.hat[i] <- as.matrix(x2optimal[i,1:length(SQEcoeff_stage2)])%*%as.matrix(SQEcoeff_stage2)
  }
  
  # create weight matrix
  SQEwt1 <- diag(length(SQEYk[[1]]))
  
  # create pseudo outcomes
  SQEPYk <- SQEYk[[1]] + stage2Q.hat
  
  # solve for betas stage 1
  x1ME <- as.matrix(SQEHk[[1]])
  x1TE <- as.matrix(SQEHk[[1]]*SQEAk[[1]])
  x1 <- cbind(x1TE, x1ME)
  SQEcoeff_stage1 <- ls.solver(x1, SQEwt1, SQEwt1, SQEPYk)  # tx effects coeff length of Hk[[1]]
  betas_stage1 <- as.matrix(SQEcoeff_stage1[1:dim(x1ME)[2]]) # stage 1 tx rule

  # create list with coefficients for both stages
  s1coeff <- as.data.frame(t(betas_stage1))
  s2coeff <- as.data.frame(t(betas_stage2))
  SQEAlgorithmResults <- list(s1coeff, s2coeff)
  names(SQEAlgorithmResults) <- c("Stage 1 Coefficients", "Stage 2 Coefficients")
  
  return(SQEAlgorithmResults)
  
}

# function to run the multi-stage augmented Q-learning estimator 
# takes in matrices for different treatment groups (mu1var, mu0var), outcome vector (Yk), treatment indicator (txind), list of history matrices for each stage (Hk), list of outcome vectors
#  at each stage (Yk), probabilities of treatment assignments (pi1k, pi0k), trial indicator (rctind), list of matrices of variables for estimating main effects at each stage (MEk),
#  list of outcomes at each stage for the trial-only (RCTy), list of treatment assignments at each stage for the trial-only (RCTa), dataset(s) with which to model potential outcome parameters 
#  (group), weight for trial data contribution to contrast term (omegaw)
# returns DTR (i.e., list of regression coefficients for each stage corresponding to the covariates from the input)

MAQE <- function(mu1var, mu0var, txind, Hk, Yk, pi1k, pi0k, rctind, MEk, RCTy, RCTa, group, omegaw){
  
  nTrial <- length(MEk[[2]][,1])  # number of rows in trial data = n
  nplusm <- length(Yk[[2]])  # number of rows in combined data = n+m
  mOS <- nplusm - nTrial  # number of rows in OS data
 
  q <- nTrial/nplusm  # proportion of subjects in trial = n/(n+m)
  minusq <- 1-q
  
  ################################################
  
  # the w1 weight matrix is what is used to square the X matrix
  # the w2 weight matrix assigns weights to outcomes
  # in this script the MAQE is programmed such that q=1 in all cases, and thus we have removed q from the description of the algorithm in the paper
  # however, at the user's discretion the weight matrices can be altered such that the outcome vectors can be weighted differently for the trial and OS
  
  w1 <- diag(nplusm) # w1 is a (n+m)x(n+m) matrix
  q1 <- rep(1/nTrial, nTrial) # 1/n for first n elements
  invm <- rep(0, mOS) # 0 for next m elements
  diag(w1) <- c(q1, invm)
  
  
  w2 <- diag(nplusm)
  qwt <- rep(1/n, nTrial) # create vector of 1/n for length of n
  minusqwt <- rep(1/mOS, mOS) # create vector of 1/m for length of m
  diag(w2) <- c(qwt, minusqwt) # w2 is (n+m)x(n+m) matrix with 1/n for first n diagonal elements and 1/m for next m diagonal elements
  
  ### stage 2  
  
  dataframe2 <- as.data.frame(cbind(Hk[[2]], Yk[[2]], txind[[2]]))
    # get mu hats stage 2
    muhats2 <- muhats.var(mu1var[[2]], mu0var[[2]], Yk[[2]], txind[[2]], rctind, group)

  # create R hats stage 2
  i=0
  Rhat2 <- rep(0,length(Yk[[2]]))
  for(i in 1:dim(dataframe2)[1]){
    if(rctind[i]==1){
      # construct R hat for trial participants
      Rhat2[i] <- ((txind[[2]][i]/pi1k[[2]][i])*(Yk[[2]][i] - muhats2[i,1])-((1-txind[[2]][i])/pi0k[[2]][i])*(Yk[[2]][i] - muhats2[i,2]) + omegaw*(muhats2[i,1] - muhats2[i,2]))
    }else if(rctind[i]==0){
      # construct R hat for OS participants
      Rhat2[i] <- (1-omegaw)*(muhats2[i,1] - muhats2[i,2])
    }
  }
  
  # solve for betas stage 2 (tx effects)
  intercept <- rep(1,length(Yk[[2]]))
  x2 <- as.matrix(as.data.frame(cbind(intercept, Hk[[2]])))
  betas_stage2 <- ls.solver(as.matrix(x2), w1, w2, Rhat2)
  
  # solve for gammas stage 2 (main effects)
  
  x2.trial <- MEk[[2]]
  x2.trialt <-t(x2.trial)
  
  
  #create an estimated value which is observed Y-A(beta_hat * xi) for each trial participant, and call this C hat
  i=0
  c.hat2 <- rep(0, nTrial)
  for(i in 1:nTrial){
    c.hat2[i] <- RCTy[[2]][i] - RCTa[[2]][i]*(x2.trial[i, 1:dim(MEk[[2]])[2]]%*%betas_stage2)
  }
  
  w.t <- diag(nTrial)
  gammas_stage2 <- ls.solver(x2.trial, w.t, w.t, c.hat2)
  
  # find optimal stage 2 tx for each subject
  
  stage2.txrule <- x2%*%betas_stage2
  stage2.optimal <- ifelse(stage2.txrule > 0, 1, 0)
  
  # calculate stage 2 Q hats based on optimal tx
  
  stage2.coeff <- c(gammas_stage2, betas_stage2)
  # then make sure to order covariates in the order that corresponds to the vector of regression coefficients
  # gammas stage 2 are main effects so
  # intercept, x1, x2, x3, resp, a1
  # betas stage 2 are interaction terms with a2 so
  # a2, x1*a2, x2*a2, x3*a2, resp*a2, a1*a2
  x2.df <- as.data.frame(x2)
  xmatrix2.trial <- cbind(x2.df, x2.df*stage2.optimal)
  
  i=0
  stage2Q.hat <- rep(0, nplusm)
  for(i in 1:nplusm){
    stage2Q.hat[i] <- as.matrix(xmatrix2.trial[i,1:length(stage2.coeff)])%*%as.matrix(stage2.coeff)
  }
  
  ### stage 1  
  
  # get mu hats stage 1
  
  # create pseudo outcome for stage 1 to pass as argument into muhats.var function
  
  psuedo.outcome <- stage2Q.hat + Yk[[1]]
  
  dataframe1 <- as.data.frame(cbind(Hk[[1]], psuedo.outcome, txind[[1]]))
    muhats1 <- muhats.var(mu1var[[1]], mu0var[[1]], psuedo.outcome, txind[[1]], rctind, group)

  # create R hats stage 1
  i=0
  Rhat1 <- rep(0,length(Yk[[1]]))
  for(i in 1:dim(dataframe1)[1]){
    if(rctind[i]==1){
      # construct R hat for trial participants
      Rhat1[i] <- ((txind[[1]][i]/pi1k[[1]][i])*(psuedo.outcome[i] - muhats1[i,1])-((1-txind[[1]][i])/pi0k[[1]][i])*(psuedo.outcome[i] - muhats1[i,2]) + omegaw*(muhats1[i,1] - muhats1[i,2]))
    }else if(rctind[i]==0){
      # construct R hat for OS participants
      Rhat1[i] <- (1-omegaw)*(muhats1[i,1] - muhats1[i,2])
    }
  }
  # solve for betas stage 1
  intercept <- rep(1,length(Yk[[1]]))
  x1 <- as.matrix(as.data.frame(cbind(intercept, Hk[[1]])))
  betas_stage1 <- ls.solver(x1, w1, w2, Rhat1)
  
  # create list with coefficients for both stages
  s1coeff <- as.data.frame(t(betas_stage1))
  s2coeff <- as.data.frame(t(betas_stage2))
  AugMethodAlgorithmResults <- list(s1coeff, s2coeff)
  names(AugMethodAlgorithmResults) <- c("Stage 1 Coefficients", "Stage 2 Coefficients")
  
  return(AugMethodAlgorithmResults)
  
  
}