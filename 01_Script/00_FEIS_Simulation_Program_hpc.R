###########################
#### Utility Functions ####
###########################


### Help function for dopar output
comb<-function(x){
  lapply(seq_along(x[[1]]), function(i) 
    do.call(Map, c(f = rbind, lapply(x, `[[`, i))))
}

comb2<-function(x){
  lapply(seq_along(x[[1]]), function(i) 
    do.call(rbind, lapply(x, `[[`, i)))
}



#############################
#### Simulation Function ####
#############################


fesim<-function(N=300, time=10, R=30, seed=123579, rob=FALSE, crs=NA,
                tol=1e-10, cov_emp=TRUE, lag=FALSE, lagprop=0.5, lagest=FALSE,
                unbal=NULL, selectivity=0,
                b1=0.3, bx1a1=1,  
                cov_a2_bx1w=0.85,
                bx1w_m=0.8, bx1w_sd=1, 
                a1_m=0, a1_sd=0, a2_m=0, a2_sd=1, 
                w_m=1, w_sd=1, x1_m=0, x1_sd=1, u_sd=0.1){
	
  
  ### Save parameters
  parameters<-c(N=N, time=time, R=R, seed=seed, rob=rob, crs=crs,
                tol=tol, cov_emp=cov_emp, lag=lag, lagprop=lagprop, lagest=lagest,
                unbal=unbal, selectivity=selectivity,
                b1=b1, bx1a1=bx1a1, 
                cov_a2_bx1w=cov_a2_bx1w, 
                bx1w_m=bx1w_m, bx1w_sd=bx1w_sd, 
                a1_m=a1_m, a1_sd=a1_sd, a2_m=a2_m, a2_sd=a2_sd, 
                w_m=w_m, w_sd=w_sd, x1_m=x1_m, x1_sd=x1_sd, u_sd=u_sd)
  
  theta=c(b1)
  
  #########################
  ### Set up basic data ###
  #########################
  
  id<-seq_along(1:N)
  id<-rep(id, each=time)
  
  t<-seq_along(1:time)
  t<-rep(t, times=N)
  
  df<-data.frame(id, t)
  
  
  
  # ##############################
  # ### Set up result matrices ###
  # ##############################
  # 
  # beta_re.df<-matrix(NA, ncol=2, nrow=R)
  # beta_fe.df<-matrix(NA, ncol=2, nrow=R)
  # beta_feis.df<-matrix(NA, ncol=2, nrow=R)
  # colnames(beta_re.df)<-colnames(beta_fe.df)<- colnames(beta_feis.df)<-c("x1", "s")
  # 
  # sd_re.df<-matrix(NA, ncol=2, nrow=R)
  # sd_fe.df<-matrix(NA, ncol=2, nrow=R)
  # sd_feis.df<-matrix(NA, ncol=2, nrow=R)
  # colnames(sd_re.df)<-colnames(sd_fe.df)<- colnames(sd_feis.df)<-c("x1", "s")
  # 
  # 
  # test.df<-matrix(NA, ncol=3, nrow=R)
  # test_p.df<-matrix(NA, ncol=3, nrow=R)
  # colnames(test.df)<- colnames(test_p.df)<-c("H0","H1", "H2")
  # 
  # creis<- vector("list", R) 
  # 
  # theta_re.df<-matrix(NA, ncol=1, nrow=R)
  # Va_re.df<-matrix(NA, ncol=1, nrow=R)
  # Ve_re.df<-matrix(NA, ncol=1, nrow=R)
  # 
  # param_emp.df <- matrix(NA, ncol=5, nrow=R)
  # colnames(param_emp.df) <- c("Var_w", "Cov_delta_a2", "Var_delta","Var_v", "biasFE_pred")
  #Var_w.df<-matrix(NA, ncol=3, nrow=R)
  #Cov_delta_a2.df<-matrix(NA, ncol=3, nrow=R)
  #Var_delta.df<-matrix(NA, ncol=3, nrow=R)
  #Var_v.df<-matrix(NA, ncol=3, nrow=R)
  #colnames(beta_re.df)<-colnames(beta_fe.df)<- colnames(beta_feis.df)<-c("x1", "x2", "s")
  
  
  ########################
  ### Start Simulation ###
  ########################
  
  set.seed(seed)
  #maxint<-.Machine$integer.max
  maxint<-2147483647
  seeds<-runif(R, min=-maxint, max=maxint)
  
  # library(doParallel)
  
  ### Register Cores
  if(is.na(crs)){
    crs<-detectCores(all.tests = FALSE, logical = TRUE)
    print(crs)
  }
  cl <- makeCluster(crs)
  registerDoParallel(cl)
  
  # for(i in 1:R){
  res <- foreach(i=1:R) %dopar% {
    
    library(Formula)
    library(matrixcalc)
    library(plm)
    library(aod)
    library(feisr)
    
    # Source feisr
    # source("00_feis.R", local = TRUE)
    # source("00_feistest.R", local = TRUE)
    # source("00_summary_functions.R", local = TRUE)
    # source("00_utility_functions.R", local = TRUE)
    
    
    set.seed(seeds[i])
    
    ### Create a slope variable
    df$w <- as.vector(replicate(N, rnorm(time, mean=w_m, sd=w_sd)))
    
    ### Create constant alpha
    df$a1 <- rep(rnorm(N, mean=a1_m, sd=a1_sd), each=time)
    
    ### Create slope alpha and coef between x and t
    M <- matrix(c(a2_sd^2, cov_a2_bx1w, 
                  cov_a2_bx1w, bx1w_sd^2), nrow=2, ncol=2)
    mu <- c(a2_m, bx1w_m)
    
    # mvrnorm for correlated norm dist
    # check whether M is positive definite
    if (is.positive.definite(M, tol=1e-8)=="TRUE") { 
      M1 <- M
      dat1 <- MASS::mvrnorm(n=N, mu=mu, Sigma=M, empirical=cov_emp)
    }
    if (is.positive.definite(M, tol=1e-8)!="TRUE") { 
      M1 <- M
      M <- Matrix::nearPD(M, conv.tol = 1e-8)
      M <- as.matrix(M$mat)
      dat1 <- MASS::mvrnorm(n=N, mu=mu, Sigma=M, empirical=cov_emp)
      #Higham, Nick (2002) Computing the nearest correlation matrix - a problem from finance; IMA Journal of Numerical Analysis 22, 329-343.
    }
    
    # Extract alpha 2
    df$a2 <- rep(dat1[,1], each=time)
    
    # Create parameter coefficient x und w
    df$bx1w <- rep(dat1[,2], each=time)
    
    ### Gen X
    df$v1 <- rnorm(N*time, mean=x1_m, sd=x1_sd)
    df$x1 <- df$v1 + bx1a1*df$a1 + df$bx1w*df$w
    
    ### Gen lag
    if(lag == TRUE){
      df <- df[order(df$id, df$w), ]
      
      df$lx1 <- ave(df$x1, by=df$id,
                    FUN = function(z) dplyr::lag(z, n=1, default=0))
      df$l2x1 <- ave(df$x1, by=df$id,
                     FUN = function(z) dplyr::lag(z, n=2, default=0))
      
    }
    
    ### Gen u
    df$u <- rnorm(N*time, mean=0, sd=u_sd) 
    
    ### Gen Y
    if(lag == TRUE){
      df$Y <- (b1*df$x1 + 
                 lagprop*b1*df$lx1 + 
                 lagprop*b1*df$l2x1 + 
                 df$a1 + 
                 df$a2*df$w +
                 df$u)
      
    }else{
      df$Y <- (b1*df$x1 + 
                 df$a1 + 
                 df$a2*df$w +
                 df$u)  
    }
    
    
    
    
    ### Implement panel selectivity
    
    if(!is.null(unbal)){
      
      if(unbal == 1){
        ### Selection 1 (on individual Y deviations)
        # Individual deviations
        df$Ydm <- ave(df$Y, df$id, FUN = function(z) z - mean(z))
        
        # Model as increasing funktion of w
        p <- selectivity
        df$z <- (p * df$Ydm / sd(df$Ydm) + (1-p) *  rnorm(N*time, mean = 0, sd = 1))
        
        df$keep <- 0
        df$keep[df$z<=0] <- 1 
        
        # set non-observed cases to NA
        df$Yorig <- df$Y
        df$Y[df$keep == 0] <- NA
      }
      
      if(unbal == 2){
        ### Selection 2 (two groups with heterogeneous effects)
        df$group <- rep(sample(c(0,1), size = N, replace = TRUE), each = time)
        prop <- mean(df$group)
        
        # Split up effect
        b1_g1 <- (2 * b1) / (2 * prop + (1 - prop))
        b1_g0 <- b1_g1 / 2
        
        df$Y <- (b1_g0 * df$x1 + 
                   (b1_g1 - b1_g0) * df$x1 * df$group +
                   df$a1 + 
                   df$a2 * df$w +
                   df$u)  
        
        ### Chance of contributing only two periods
        p <- selectivity
        df$pr <- 1 - p
        df$pr[df$group == 0] <- p
        
        df$twoper <- rep(rbinom(N, 1, prob = df$pr[df$t == 1]), each = time)
        
        # select two time point randomly in selected groups/individuals
        sample2per <- function(z){
          z[] <- 0
          z[sample(1:length(z), 2)] <- 1
          return(z)
        }
        
        df$keep <- 1
        df$keep[df$twoper == 1] <- ave(df$keep[df$twoper == 1], df$id[df$twoper == 1], 
                                       FUN = function(z) sample2per(z)) 
        
        # set non-observed cases to NA
        df$Yorig <- df$Y
        df$Y[df$keep == 0] <- NA
      }
    }
    
    
    
    df$t2<-df$t # For including as continuous in plm
    
    
    # ### Store empirical parameters (to calculate predicted FE bias)   
    # 
    # param_emp.df[i,1] <- var(df$w-ave(df$w, by=df$id, FUN=mean))
    # param_emp.df[i,2] <- cov((df$bx1w-mean(df$bx1w)),(df$a2-mean(df$a2)))
    # #param_emp.df[i,2] <- cov(df$bx1w,(df$a2-mean(df$a2)))
    # param_emp.df[i,3] <- var(unique(df$bx1w-mean(df$bx1w)))
    # param_emp.df[i,4] <- var(df$v1-ave(df$v1, by=df$id, FUN=mean))
    # param_emp.df[i,5] <- ((param_emp.df[i,1] * param_emp.df[i,2]) / 
    #                         ((param_emp.df[i,1] * param_emp.df[i,3]) + param_emp.df[i,4]))
    
    ################
    ### Analysis ###
    ################
    
    if(lagest == FALSE){
      ### RE
      re.mod <- plm(Y ~ x1 + w, data=df, index=c("id", "t"),
                    effect="individual", model="random")
      
      ### FE
      fe.mod <- plm(Y ~ x1 + w, data=df, index=c("id", "t"),
                    effect="individual", model="within")
      
      ### FEIS
      if(!is.null(unbal)){
        suppressWarnings( # Suppress "drop groups" warning if unbalanced
          feis.mod<-feis("Y ~ x1  | w", data=df, id="id")
        )
      }else{
        feis.mod<-feis("Y ~ x1  | w", data=df, id="id")
      }
    }
    
    if(lagest == TRUE){
      ### RE
      re.mod <- plm(Y ~ x1 + lx1 + l2x1 + w, data=df, index=c("id", "t"),
                    effect="individual", model="random")
      
      ### FE
      fe.mod <- plm(Y ~ x1 + lx1 + l2x1 + w, data=df, index=c("id", "t"),
                    effect="individual", model="within")
      
      ### FEIS
      if(!is.null(unbal)){
        suppressWarnings( # Suppress "drop groups" warning if unbalanced
          feis.mod<-feis("Y ~ x1 + lx1 + l2x1 | w", data=df, id="id")
        )
      }else{
        feis.mod<-feis("Y ~ x1 + lx1 + l2x1 | w", data=df, id="id")
      }
    }
    
    
    ### Augmented Regression Test
    hf<-feistest(feis.mod, type="all", robust=rob)
    
    
    ### Set up results mat
    if(lagest == FALSE){
      beta_re.df<-matrix(NA, ncol=2, nrow=1)
      beta_fe.df<-matrix(NA, ncol=2, nrow=1)
      beta_feis.df<-matrix(NA, ncol=2, nrow=1)
      colnames(beta_re.df)<-colnames(beta_fe.df)<- colnames(beta_feis.df)<-c("x1", "s")
      
      sd_re.df<-matrix(NA, ncol=2, nrow=1)
      sd_fe.df<-matrix(NA, ncol=2, nrow=1)
      sd_feis.df<-matrix(NA, ncol=2, nrow=1)
      colnames(sd_re.df)<-colnames(sd_fe.df)<- colnames(sd_feis.df)<-c("x1", "s")
    }
    
    if(lagest == TRUE){
      beta_re.df <- matrix(NA, ncol=4, nrow=1)
      beta_fe.df <- matrix(NA, ncol=4, nrow=1)
      beta_feis.df <- matrix(NA, ncol=4, nrow=1)
      colnames(beta_re.df) <- colnames(beta_fe.df) <-  colnames(beta_feis.df) <- c("x1", "lx1", "l2x1", "s")
      
      sd_re.df <- matrix(NA, ncol=4, nrow=1)
      sd_fe.df <- matrix(NA, ncol=4, nrow=1)
      sd_feis.df <- matrix(NA, ncol=4, nrow=1)
      colnames(sd_re.df) <- colnames(sd_fe.df) <-  colnames(sd_feis.df) <- c("x1", "lx1", "l2x1", "s")
    }
    
    
    test.df<-matrix(NA, ncol=3, nrow=1)
    test_p.df<-matrix(NA, ncol=3, nrow=1)
    colnames(test.df)<- colnames(test_p.df)<-c("H0","H1", "H2")
    
    creis<- vector("list", 1) 
    
    theta_re.df<-matrix(NA, ncol=1, nrow=1)
    Va_re.df<-matrix(NA, ncol=1, nrow=1)
    Ve_re.df<-matrix(NA, ncol=1, nrow=1)
    
    param_emp.df <- matrix(NA, ncol=5, nrow=1)
    colnames(param_emp.df) <- c("Var_w", "Cov_delta_a2", "Var_delta","Var_v", "biasFE_pred")
    
    
    ### Paste values
    beta_re.df[1,]<-re.mod$coefficients[-1]
    beta_fe.df[1,]<-fe.mod$coefficients
    beta_feis.df[1,1:(ncol(beta_feis.df)-1)]<-feis.mod$coefficients
    
    sd_re.df[1,]<-sqrt(diag(vcov(re.mod)))[-1]
    sd_fe.df[1,]<-sqrt(diag(vcov(fe.mod)))
    sd_feis.df[1,1:(ncol(sd_feis.df)-1)]<-sqrt(diag(feis.mod$vcov))
    
    test.df[1,1]<-hf$wald_feis$result$chi2[1]
    test_p.df[1,1]<-hf$wald_feis$result$chi2[3]
    
    test.df[1,2]<-hf$wald_fe$result$chi2[1]
    test_p.df[1,2]<-hf$wald_fe$result$chi2[3]
    
    test.df[1,3]<-hf$wald_re$result$chi2[1]
    test_p.df[1,3]<-hf$wald_re$result$chi2[3]
    
    creis[[1]]<-hf$CREIS
    
    theta_re.df[1,] <- mean(re.mod$ercomp$theta)
    Va_re.df[1,] <- re.mod$ercomp$sigma2[1]
    Ve_re.df[1,] <- re.mod$ercomp$sigma2[2]
    
    param_emp.df[1,1] <- var(df$w-ave(df$w, by=df$id, FUN=mean))
    param_emp.df[1,2] <- cov((df$bx1w-mean(df$bx1w)),(df$a2-mean(df$a2)))
    # param_emp.df[1,2] <- cov(df$bx1w,(df$a2-mean(df$a2)))
    param_emp.df[1,3] <- var(unique(df$bx1w-mean(df$bx1w)))
    param_emp.df[1,4] <- var(df$v1-ave(df$v1, by=df$id, FUN=mean))
    param_emp.df[1,5] <- ((param_emp.df[1,1] * param_emp.df[1,2]) / 
                            ((param_emp.df[1,1] * param_emp.df[1,3]) + param_emp.df[1,4]))
    
    
    
    # print(seeds[i])
    if(i==1){
      cat(paste("Simulations completed (by 10):"," "))
    } 
    if(i %% 50==0){
      cat(paste(" ", i," "))
    }else if(i %% 10==0){
      cat(paste("+"))
    }
    if(i==R){
      cat("\n")
    }
    
    res <- list(beta_re.df      = beta_re.df,
                beta_fe.df      = beta_fe.df,
                beta_feis.df    = beta_feis.df,
                sd_re.df      = sd_re.df,
                sd_fe.df      = sd_fe.df,
                sd_feis.df    = sd_feis.df,
                test.df		  = test.df,
                test_p.df	  = test_p.df,
                #creis         = creis,
                param_emp.df  = param_emp.df,
                theta_re.df   = theta_re.df,
                Va_re.df      = Va_re.df,
                Ve_re.df      = Ve_re.df)
    
    return(res)
    
  }
  
  ### Stop Cluster
  
  stopCluster(cl)
  
  ### Reshape results
  
  res <- comb2(res)
  
  ### Create output element
  
  result<-list(parameters    = parameters,
               theta         = theta,
               lagprop       = lagprop,
               seeds         = seeds,
               beta_re.df    = res[[1]],
               beta_fe.df    = res[[2]],
               beta_feis.df  = res[[3]],
               sd_re.df      = res[[4]],
               sd_fe.df      = res[[5]],
               sd_feis.df    = res[[6]],
               test.df       = res[[7]],
               test_p.df     = res[[8]],
               #CREIS         = res[[9]],
               param_emp.df  = res[[9]],
               theta_re.df   = res[[10]],
               Va_re.df      = res[[11]],
               Ve_re.df      = res[[12]])
  
  class(result) <- c("fesim")
  
  return(result)
}



#######################################
#### Simulation Function Bootstrap ####
#######################################


fesim_bs<-function(N=300, time=10, R=30, seed=123579, bsR=10, crs=NA,
                   tol=1e-10, cov_emp=TRUE,
                   b1=0.3, bx1a1=1,  
                   cov_a2_bx1w=0.85,
                   bx1w_m=0.8, bx1w_sd=1, 
                   a1_m=0, a1_sd=0, a2_m=0, a2_sd=1, 
                   w_m=1, w_sd=1, x1_m=0, x1_sd=1, u_sd=0.1){
  
  
  ### Save parameters
  parameters<-c(N=N, time=time, R=R, seed=seed, bsR=10, crs=NA,
                tol=tol, cov_emp=cov_emp,
                b1=b1, bx1a1=bx1a1, 
                cov_a2_bx1w=cov_a2_bx1w, 
                bx1w_m=bx1w_m, bx1w_sd=bx1w_sd, 
                a1_m=a1_m, a1_sd=a1_sd, a2_m=a2_m, a2_sd=a2_sd, 
                w_m=w_m, w_sd=w_sd, x1_m=x1_m, x1_sd=x1_sd, u_sd=u_sd)
  
  theta=c(b1)
  
  #########################
  ### Set up basic data ###
  #########################
  
  id<-seq_along(1:N)
  id<-rep(id, each=time)
  
  t<-seq_along(1:time)
  t<-rep(t, times=N)
  
  df<-data.frame(id, t)
  
  
  
  # ##############################
  # ### Set up result matrices ###
  # ##############################
  # 
  # beta_re.df<-matrix(NA, ncol=2, nrow=R)
  # beta_fe.df<-matrix(NA, ncol=2, nrow=R)
  # beta_feis.df<-matrix(NA, ncol=2, nrow=R)
  # colnames(beta_re.df)<-colnames(beta_fe.df)<- colnames(beta_feis.df)<-c("x1", "s")
  # 
  # sd_re.df<-matrix(NA, ncol=2, nrow=R)
  # sd_fe.df<-matrix(NA, ncol=2, nrow=R)
  # sd_feis.df<-matrix(NA, ncol=2, nrow=R)
  # colnames(sd_re.df)<-colnames(sd_fe.df)<- colnames(sd_feis.df)<-c("x1", "s")
  # 
  # 
  # test.df<-matrix(NA, ncol=3, nrow=R)
  # test_p.df<-matrix(NA, ncol=3, nrow=R)
  # colnames(test.df)<- colnames(test_p.df)<-c("H0","H1", "H2")
  # 
  # creis<- vector("list", R) 
  # 
  # theta_re.df<-matrix(NA, ncol=1, nrow=R)
  # Va_re.df<-matrix(NA, ncol=1, nrow=R)
  # Ve_re.df<-matrix(NA, ncol=1, nrow=R)
  # 
  # param_emp.df <- matrix(NA, ncol=5, nrow=R)
  # colnames(param_emp.df) <- c("Var_w", "Cov_delta_a2", "Var_delta","Var_v", "biasFE_pred")
  #Var_w.df<-matrix(NA, ncol=3, nrow=R)
  #Cov_delta_a2.df<-matrix(NA, ncol=3, nrow=R)
  #Var_delta.df<-matrix(NA, ncol=3, nrow=R)
  #Var_v.df<-matrix(NA, ncol=3, nrow=R)
  #colnames(beta_re.df)<-colnames(beta_fe.df)<- colnames(beta_feis.df)<-c("x1", "x2", "s")
  
  
  ########################
  ### Start Simulation ###
  ########################
  
  set.seed(seed)
  #maxint<-.Machine$integer.max
  maxint<-2147483647
  seeds<-runif(R, min=-maxint, max=maxint)
  
  # library(doParallel)
  
  ### Register Cores
  if(is.na(crs)){
    crs<-detectCores(all.tests = FALSE, logical = TRUE)
    print(crs)
  }
  cl <- makeCluster(crs)
  registerDoParallel(cl)
  
  # for(i in 1:R){
  res <- foreach(i=1:R) %dopar% {
    
    library(Formula)
    library(matrixcalc)
    library(plm)
    library(aod)
    library(feisr)
    
    # Source feisr
    # source("00_feis.R", local = TRUE)
    # source("00_feistest.R", local = TRUE)
    # source("00_summary_functions.R", local = TRUE)
    # source("00_utility_functions.R", local = TRUE)
    
    
    set.seed(seeds[i])
    
    ### Create a slope variable
    df$w <- as.vector(replicate(N, rnorm(time, mean=w_m, sd=w_sd)))
    
    ### Create constant alpha
    df$a1 <- rep(rnorm(N, mean=a1_m, sd=a1_sd), each=time)
    
    ### Create slope alpha and coef between x and t
    M <- matrix(c(a2_sd^2, cov_a2_bx1w, 
                  cov_a2_bx1w, bx1w_sd^2), nrow=2, ncol=2)
    mu <- c(a2_m, bx1w_m)
    
    # mvrnorm for correlated norm dist
    # check whether M is positive definite
    if (is.positive.definite(M, tol=1e-8)=="TRUE") { 
      M1 <- M
      dat1 <- MASS::mvrnorm(n=N, mu=mu, Sigma=M, empirical=cov_emp)
    }
    if (is.positive.definite(M, tol=1e-8)!="TRUE") { 
      M1 <- M
      M <- Matrix::nearPD(M, conv.tol = 1e-8)
      M <- as.matrix(M$mat)
      dat1 <- MASS::mvrnorm(n=N, mu=mu, Sigma=M, empirical=cov_emp)
      #Higham, Nick (2002) Computing the nearest correlation matrix - a problem from finance; IMA Journal of Numerical Analysis 22, 329-343.
    }
    
    # Extract alpha 2
    df$a2 <- rep(dat1[,1], each=time)
    
    # Create parameter coefficient x und w
    df$bx1w <- rep(dat1[,2], each=time)
    
    ### Gen X
    df$v1 <- rnorm(N*time, mean=x1_m, sd=x1_sd)
    df$x1 <- df$v1 + bx1a1*df$a1 + df$bx1w*df$w
    
    ### Gen u
    df$u <- rnorm(N*time, mean=0, sd=u_sd) 
    
    ### Gen Y
    df$Y <- (b1*df$x1 + 
               df$a1 + 
               df$a2*df$w +
               df$u)
    
    
    
    df$t2<-df$t # For including as continuous in plm
    
    
    ### Store empirical parameters (to calculate predicted FE bias)   
    
    # param_emp.df[i,1] <- var(df$w-ave(df$w, by=df$id, FUN=mean))
    # param_emp.df[i,2] <- cov((df$bx1w-mean(df$bx1w)),(df$a2-mean(df$a2)))
    # #param_emp.df[i,2] <- cov(df$bx1w,(df$a2-mean(df$a2)))
    # param_emp.df[i,3] <- var(unique(df$bx1w-mean(df$bx1w)))
    # param_emp.df[i,4] <- var(df$v1-ave(df$v1, by=df$id, FUN=mean))
    # param_emp.df[i,5] <- ((param_emp.df[i,1] * param_emp.df[i,2]) / 
    #       ((param_emp.df[i,1] * param_emp.df[i,3]) + param_emp.df[i,4]))
    
    # df$a2_dm<-df$a2-mean(df$a2)
    # df$w_dm<-df$w-ave(df$w, by=df$id, FUN=function(x) mean(x))
    # df$x1_dm<-df$x1-ave(df$x1, by=df$id, FUN=function(x) mean(x))
    # df$v1_dm<-df$v1-ave(df$v1, by=df$id, FUN=function(x) mean(x))
    # df$Y_dm<-df$Y-ave(df$Y, by=df$id, FUN=function(x) mean(x))
    # df$u_dm<-df$u-ave(df$u, by=df$id, FUN=function(x) mean(x))
    # 
    # df$bx1w_dm<-df$bx1w-mean(df$bx1w)
    # 
    # 
    # # bias 
    # (var(df$w_dm)*cov(df$a2_dm, df$bx1w)) / (var(df$w_dm)*(var(unique(df$bx1w_dm))) + var(df$v1_dm))
    #     
    
    ################
    ### Analysis ###
    ################
    
    ### RE
    re.mod<-plm(Y ~ x1 + w, data=df, index=c("id", "t"),
                effect="individual", model="random")
    
    ### FE
    fe.mod<-plm(Y ~ x1 + w, data=df, index=c("id", "t"),
                effect="individual", model="within")
    
    ### FEIS
    feis.mod<-feis("Y ~ x1  | w", data=df, id="id")
    
    ### Augmented Regression Test
    hf<-bsfeistest(feis.mod, type="all", rep=bsR, prog = F)
    
    
    ### Set up results mat
    beta_re.df<-matrix(NA, ncol=2, nrow=1)
    beta_fe.df<-matrix(NA, ncol=2, nrow=1)
    beta_feis.df<-matrix(NA, ncol=2, nrow=1)
    colnames(beta_re.df)<-colnames(beta_fe.df)<- colnames(beta_feis.df)<-c("x1", "s")
    
    sd_re.df<-matrix(NA, ncol=2, nrow=1)
    sd_fe.df<-matrix(NA, ncol=2, nrow=1)
    sd_feis.df<-matrix(NA, ncol=2, nrow=1)
    colnames(sd_re.df)<-colnames(sd_fe.df)<- colnames(sd_feis.df)<-c("x1", "s")
    
    
    test.df<-matrix(NA, ncol=3, nrow=1)
    test_p.df<-matrix(NA, ncol=3, nrow=1)
    colnames(test.df)<- colnames(test_p.df)<-c("H0","H1", "H2")
    
    creis<- vector("list", 1) 
    
    theta_re.df<-matrix(NA, ncol=1, nrow=1)
    Va_re.df<-matrix(NA, ncol=1, nrow=1)
    Ve_re.df<-matrix(NA, ncol=1, nrow=1)
    
    param_emp.df <- matrix(NA, ncol=5, nrow=1)
    colnames(param_emp.df) <- c("Var_w", "Cov_delta_a2", "Var_delta","Var_v", "biasFE_pred")
    
    
    ### Paste values
    beta_re.df[1,]<-re.mod$coefficients[-1]
    beta_fe.df[1,]<-fe.mod$coefficients
    beta_feis.df[1,1]<-feis.mod$coefficients
    
    sd_re.df[1,]<-sqrt(diag(vcov(re.mod)))[-1]
    sd_fe.df[1,]<-sqrt(diag(vcov(fe.mod)))
    sd_feis.df[1,1]<-sqrt(diag(feis.mod$vcov))
    
    test.df[1,1]<-hf$wald_feis$result$chi2[1]
    test_p.df[1,1]<-hf$wald_feis$result$chi2[3]
    
    test.df[1,2]<-hf$wald_fe$result$chi2[1]
    test_p.df[1,2]<-hf$wald_fe$result$chi2[3]
    
    test.df[1,3]<-hf$wald_re$result$chi2[1]
    test_p.df[1,3]<-hf$wald_re$result$chi2[3]
    
    creis[[1]]<-hf$CREIS
    
    theta_re.df[1,] <- mean(re.mod$ercomp$theta)
    Va_re.df[1,] <- re.mod$ercomp$sigma2[1]
    Ve_re.df[1,] <- re.mod$ercomp$sigma2[2]
    
    param_emp.df[1,1] <- var(df$w-ave(df$w, by=df$id, FUN=mean))
    param_emp.df[1,2] <- cov((df$bx1w-mean(df$bx1w)),(df$a2-mean(df$a2)))
    # param_emp.df[1,2] <- cov(df$bx1w,(df$a2-mean(df$a2)))
    param_emp.df[1,3] <- var(unique(df$bx1w-mean(df$bx1w)))
    param_emp.df[1,4] <- var(df$v1-ave(df$v1, by=df$id, FUN=mean))
    param_emp.df[1,5] <- ((param_emp.df[1,1] * param_emp.df[1,2]) / 
                            ((param_emp.df[1,1] * param_emp.df[1,3]) + param_emp.df[1,4]))
    
    
    
    # print(seeds[i])
    if(i==1){
      cat(paste("Simulations completed (by 10):"," "))
    } 
    if(i %% 50==0){
      cat(paste(" ", i," "))
    }else if(i %% 10==0){
      cat(paste("+"))
    }
    if(i==R){
      cat("\n")
    }
    
    res <- list(beta_re.df      = beta_re.df,
                beta_fe.df      = beta_fe.df,
                beta_feis.df    = beta_feis.df,
                sd_re.df      = sd_re.df,
                sd_fe.df      = sd_fe.df,
                sd_feis.df    = sd_feis.df,
                test.df		  = test.df,
                test_p.df	  = test_p.df,
                #creis         = creis,
                param_emp.df  = param_emp.df,
                theta_re.df   = theta_re.df,
                Va_re.df      = Va_re.df,
                Ve_re.df      = Ve_re.df)
    
    return(res)
    
  }
  
  ### Stop Cluster
  
  stopCluster(cl)
  
  ### Reshape results
  
  res <- comb2(res)
  
  ### Create output element
  
  result<-list(parameters    = parameters,
               theta         = theta,
               seeds         = seeds,
               beta_re.df    = res[[1]],
               beta_fe.df    = res[[2]],
               beta_feis.df  = res[[3]],
               sd_re.df      = res[[4]],
               sd_fe.df      = res[[5]],
               sd_feis.df    = res[[6]],
               test.df       = res[[7]],
               test_p.df     = res[[8]],
               #CREIS         = res[[9]],
               param_emp.df  = res[[9]],
               theta_re.df   = res[[10]],
               Va_re.df      = res[[11]],
               Ve_re.df      = res[[12]])
  
  class(result) <- c("fesim")
  
  return(result)
}




####################################
#### Summary functions for Sims ####
####################################

### RMSE
rmse<-function(sim=NULL, true=NULL){
  se<-(sim-true)^2
  mse<-mean(se)
  rmse<-sqrt(mse)
  return(rmse)
}


### Summary
summary.fesim<-function(sim=NULL){
  
  ### Set up results matrix
  mat<-matrix(nrow=9, ncol=2)
  colnames(mat)<-c("Bias", "RMSE")
  rownames(mat)<-c("RE", "x1", "x2",
                   "FE", "x1", "x2",
                   "FEIS", "x1", "x2")
  
  comp<-matrix(nrow=nrow(sim$beta_fe.df), ncol=2)
  
  mat_art<-matrix(nrow=1, ncol=3)
  colnames(mat_art)<-c("H0","H1", "H2")
  
  ### Get true parameters
  b1<-as.numeric(sim$theta[1])
  b2<-as.numeric(sim$theta[2])
  
  
  ### Errors
  
  # RE
  mat[2,1]<-mean(sim$beta_re.df[,1])-b1
  mat[2,2]<-rmse(sim$beta_re.df[,1], b1)
  mat[3,1]<-mean(sim$beta_re.df[,2])-b2
  mat[3,2]<-rmse(sim$beta_re.df[,2], b2)
  
  # FE
  mat[5,1]<-mean(sim$beta_fe.df[,1])-b1
  mat[5,2]<-rmse(sim$beta_fe.df[,1], b1)
  mat[6,1]<-mean(sim$beta_fe.df[,2])-b2
  mat[6,2]<-rmse(sim$beta_fe.df[,2], b2)
  
  # FEIS
  mat[8,1]<-mean(sim$beta_feis.df[,1])-b1
  mat[8,2]<-rmse(sim$beta_feis.df[,1], b1)
  mat[9,1]<-mean(sim$beta_feis.df[,2])-b2
  mat[9,2]<-rmse(sim$beta_feis.df[,2], b2)
  
  ### Compare FEIS-FE
  comp1<-abs(sim$beta_feis.df[,1]-b1) - abs(sim$beta_fe.df[,1]-b1)
  comp2<-abs(sim$beta_feis.df[,2]-b2) - abs(sim$beta_fe.df[,2]-b2)
  
  comp[,1]<-comp1
  comp[,2]<-comp2
  
  ### Augmented Regression test
  test_p.df<-sim$test_p.df
  h1<-length(which(test_p.df[,1]<0.05))/nrow(test_p.df)
  h2<-length(which(test_p.df[,2]<0.05))/nrow(test_p.df)
  h0<-length(which(test_p.df[,3]<0.05))/nrow(test_p.df)
  
  mat_art[1,1]<-h1
  mat_art[1,2]<-h2
  mat_art[1,3]<-h0
  
  
  
  res<-list(mat=mat, mat_art=mat_art, comp=comp,
            theta=sim$theta, parameters=sim$parameters)
  
  class(res) <- c("summary.fesim")
  
  res
}


### Print Summary
print.summary.fesim <- function(x, digits = 5, width=getOption("width"), ...){
  
  
  ## Make printable mat direct
  mat<-x$mat
  mat_art<-x$mat_art
  comp<-x$comp
  
  feis_lower<-c(length(which(comp[,1]<0))/nrow(comp), length(which(comp[,2]<0))/nrow(comp))
  
  mat<-round(mat, digits)
  
  cat("\n")
  cat("\nCall:\n")
  print(x$parameters)
  cat("\n")
  cat("\nBeta:\n")
  print(x$theta)
  cat("\n")
  
  cat("\nResults Error\n")
  print(mat, row.names = FALSE)
  cat("\nFE rejected\n")
  cat(round(mat_art[1,1]*100,2), "%")
  cat("\nRE rejected (vs. FE)\n")
  cat(round(mat_art[1,2]*100,2), "%")
  cat("\nRE rejected (vs. FEIS)\n")
  cat(round(mat_art[1,3]*100,2), "%")
  cat("\nFEIS lower bias than FE\n")
  cat("x1:", round(feis_lower[1]*100,2), "%", "x2:", round(feis_lower[2]*100,2), "%")
  cat("\n")
  
  invisible(x)
}


