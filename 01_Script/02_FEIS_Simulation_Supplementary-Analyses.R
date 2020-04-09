#### ---------------------------------------------------------------####
# Replication supplementary simulation results
# R Script for article 
# Rüttenauer T. and Ludwig V. 2020. Fixed Effects Individual Slopes: 
# Accounting and Testing for Heterogeneous Effects in Panel Data 
# or other Multilevel Models. Sociological Methods and Research
#### ---------------------------------------------------------------####




##################
#### Packages ####
##################
# The following packages are needed for the simulations and the outputs:

### Install
# install.packages("Formula")
# install.packages("Matrix")
# install.packages("MASS")
# install.packages("matrixcalc")
# install.packages("plm")
# install.packages("aod")
# install.packages("feisr") 
# install.packages("doParallel")
# 
# install.packages("ggplot2")
# install.packages("lattice")
# install.packages("gridExtra")
# install.packages("ggpubr")
# install.packages("extrafont")

### Load
library(Formula)
library(Matrix)
library(MASS)
library(matrixcalc)
library(plm)
library(feisr)
library(aod)
library(doParallel)

library(ggplot2)
library(lattice)
library(gridExtra)
library(ggpubr)
library(extrafont)
loadfonts()



###################################
#### Specify working directory ####
###################################

# We recommend using one directory for R Scripts, Data, and Output
# as setting different directories may be more complicated depending 
# on the HPC operating system

# On HPC you may need to specify externally. In this case comment out setwd

setwd("")


#####################################
#### Load functions from program ####
#####################################

### Version without DoPar
source("../01_Script/00_FEIS_Simulation_Program.R", local=T)

### Version with DoPar (for hpc)
# source("../01_Script/00_FEIS_Simulation_Program_hpc.R", local=T)





####----------------------------------------------------------------------####
####                          Small N, large T                            ####
####----------------------------------------------------------------------####



############################
#### Specify parameters ####
############################

### Number of obs per sim
nobs <- 20

### Time periods per sim
tobs <- 30

### Replications per sim
reps <- 1000

### Replications per bootstrap run
bsreps <- 100

### Number of CPU cores used (for HPC verion of program functions)
# If set to NA in combination with hpc program verion, all cores are used (only recommended for hpc)
crs <- NA





#################################################
#### 1) Simulations varying phi, alpha1: ART ####
#################################################


### Varying phi1:
# covariance a2 and phi1 (parameter: cov_a2_bx1w)
phi1 <- c(-0.8,seq(-0.5,-0.3,0.1),seq(-0.25,0.25,0.05),seq(0.3,0.5,0.1),0.8)

### Varying bias due to alpha1
balpha1 <- c(0,1)

### results matrix
df0 <- data.frame(phi1=phi1)
df1 <- df0

start_time <- Sys.time()
### Simulations
for (ba1 in balpha1) {
  i=1
  for (p1 in phi1) {
    sim <- fesim(N=nobs, time=tobs, R=reps, rob=T, crs=crs, tol=1e-10, seed=1565761,
                 b1=1, bx1a1=ba1,
                 cov_a2_bx1w=p1,
                 bx1w_m=1, bx1w_sd=(1^0.5),
                 a1_m=1, a1_sd=(4^0.5), a2_m=0.1, a2_sd=(4^0.5),
                 w_m=1, w_sd=1, x1_m=1, x1_sd=(1^0.5), u_sd=1)
    if (ba1==0) {
      #relative bias
      df0$b1[i] <- mean(sim$beta_re.df[,1])
      df0$b2[i] <- mean(sim$beta_fe.df[,1])
      df0$b3[i] <- mean(sim$beta_feis.df[,1])
      df0$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df0$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df0$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df0$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df0$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df0$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df0$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df0$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df0$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df0$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      df0$theta_re[i] <- mean(sim$theta_re.df)
      df0$Va_re[i] <- mean(sim$Va_re.df)
      df0$Ve_re[i] <- mean(sim$Ve_re.df)
      df0$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    else if (ba1==1) {
      #relative bias
      df1$b1[i] <- mean(sim$beta_re.df[,1])
      df1$b2[i] <- mean(sim$beta_fe.df[,1])
      df1$b3[i] <- mean(sim$beta_feis.df[,1])
      df1$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df1$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df1$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df1$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df1$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df1$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df1$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df1$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df1$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df1$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      df1$theta_re[i] <- mean(sim$theta_re.df)
      df1$Va_re[i] <- mean(sim$Va_re.df)
      df1$Ve_re[i] <- mean(sim$Ve_re.df)
      df1$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    
    save(sim,
         file=paste0("Simulation1_", ba1, "_", i, "_smallN.RData", sep=""))
    
    i=i+1
  }
}
end_time <- Sys.time()
end_time - start_time

save(df0, df1, file="ART_smallN.RData")



#################################################
#### 2) Simulations varying phi, alpha1: BST ####
#################################################


### Varying phi1:
# covariance a2 and phi1 (parameter: cov_a2_bx1w)
phi1 <- c(-0.8,seq(-0.5,-0.3,0.1),seq(-0.25,0.25,0.05),seq(0.3,0.5,0.1),0.8)

### Varying bias due to alpha1
balpha1 <- c(0,1)

### results matrix
df0bs <- data.frame(phi1=phi1)
df1bs <- df0bs

### Simulations
for (ba1 in balpha1) {
  i=1
  for (p1 in phi1) {
    sim <- fesim(N=nobs, time=tobs, R=reps, rob=T, crs=crs, tol=1e-10, seed=1565761,
                 b1=1, bx1a1=ba1,
                 cov_a2_bx1w=p1,
                 bx1w_m=1, bx1w_sd=(1^0.5),
                 a1_m=1, a1_sd=(4^0.5), a2_m=0.1, a2_sd=(4^0.5),
                 w_m=1, w_sd=1, x1_m=1, x1_sd=(1^0.5), u_sd=1)
    if (ba1==0) {
      #relative bias
      df0bs$b1[i] <- mean(sim$beta_re.df[,1])
      df0bs$b2[i] <- mean(sim$beta_fe.df[,1])
      df0bs$b3[i] <- mean(sim$beta_feis.df[,1])
      df0bs$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df0bs$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df0bs$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df0bs$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0bs$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0bs$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0bs$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df0bs$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df0bs$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df0bs$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df0bs$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df0bs$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df0bs$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
    }
    else if (ba1==1) {
      #relative bias
      df1bs$b1[i] <- mean(sim$beta_re.df[,1])
      df1bs$b2[i] <- mean(sim$beta_fe.df[,1])
      df1bs$b3[i] <- mean(sim$beta_feis.df[,1])
      df1bs$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df1bs$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df1bs$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df1bs$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1bs$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1bs$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1bs$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df1bs$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df1bs$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df1bs$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df1bs$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df1bs$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df1bs$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
    }
    
    save(sim,
         file=paste0("Simulation1_bs_", ba1, "_", i, "_smallN.RData", sep=""))
    
    i=i+1
  }
}

save(df0bs,df1bs,file="BSHT_smallN.RData")




################################################
### Plot results 1) and 2)                  ####
################################################


# load("ART_smallN.RData")
# load("BSHT_smallN.RData")



### 1) results for alpha=0
res <- reshape(df0, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 
res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))


simplot1_1 <- ggplot(res, aes(phi1,bias1, group=time, color=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  #  geom_hline(yintercept = 0, color="red") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias of ",beta,"  +/- 2 SD"))) +
  ggtitle("") +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.5,0.6), breaks=c(seq(-0.5,-0.1,0.1), 0, seq(0.1,0.6,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.9,0.1), legend.justification=c(0.9,0.1),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot1_1


### ART Rejection rate
res <- reshape(df0, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot1_2 <- ggplot(res, aes(phi1,rej1, group=time, color=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle("") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot1_2



### BSHT Rejection rate
res <- reshape(df0bs, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
#test_p.df[i,1]<-hf$wald_feis$result$chi2[3]
#test_p.df[i,2]<-hf$wald_fe$result$chi2[3]
#test_p.df[i,3]<-hf$wald_re$result$chi2[3]
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot1_2bs <- ggplot(res, aes(phi1,rej1, group=time, color=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle("") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot1_2bs



### 2) results for alpha>0

res <- reshape(df1, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 
res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))


simplot1_3 <- ggplot(res, aes(phi1,bias1, group=time, color=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  #  geom_hline(yintercept = 0, color="red") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias of ",beta,"  +/- 2 SD"))) +
  ggtitle(" ") +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.5,0.6), breaks=c(seq(-0.5,-0.1,0.1), 0, seq(0.1,0.6,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.9,0.1), legend.justification=c(0.9,0.1),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot1_3



### ART Rejection rate
res <- reshape(df1, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot1_4 <- ggplot(res, aes(phi1,rej1, group=time, color=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle("") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot1_4



### BSHT Rejection rate
res <- reshape(df1bs, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot1_4bs <- ggplot(res, aes(phi1,rej1, group=time, color=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")")))+
  ylab("Proportion rejected") +
  ggtitle("") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot1_4bs


### Combine Figures
library(ggpubr)
Simplot1 <- ggarrange(simplot1_1,simplot1_2,simplot1_2bs,simplot1_3,simplot1_4,simplot1_4bs,
                      ncol=3, nrow=2, common.legend = F, 
                      labels=c("A) ", "B) ", "C) ", "D) ", "E) ", "F) "), font.label = list(size = 20)) 


cairo_pdf(file=paste("../03_Output/Output_", "Simplot1_smallN.pdf", sep=""), width=17.54, height=10, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
Simplot1
dev.off()








#################################################################
#### 3) Simulations varying phi, sigma1, sigma2, sigma3: ART ####
#################################################################


# Vary phi1 (Cov alpha2, delta) and sigma1 (Var delta) simultaneously : 
phi1 <- c(seq(0,0.2,0.05),seq(0.3,0.5,0.1),0.65,0.8)
sigma1 <- c(0.2,0.5,1,3,5)

# Vary sigma2 (Var of w) and sigma3 (Var of x)
sigma2 <- c(0.2,1,5)
sigma3 <- c(0.2,1,5)



bl <- list()
biasl <- list()
biassdl <- list()
rejl <- list()
rmsel <- list()
bias_predl <- list()
for (s2 in sigma2) {
  for (s3 in sigma3) {
    for (x in c(1,2,3)) {
      bla <- paste("b",x,s2,s3,sep="_") 
      bl[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
      bla <- paste("bias",x,s2,s3,sep="_") 
      biasl[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
      bla <- paste("biassd",x,s2,s3,sep="_") 
      biassdl[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
      bla <- paste("rej",x,s2,s3,sep="_") 
      rejl[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
      bla <- paste("rmse",x,s2,s3,sep="_") 
      rmsel[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
      bla <- paste("bias_pred",x,s2,s3,sep="_") 
      bias_predl[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
    }
  }    
}  



### Run simulations

i=1
for (p1 in phi1) {
  j=1
  for (s1 in sigma1) {
    k=1
    for (s2 in sigma2) {
      l=1
      for (s3 in sigma3) {
        
        sim <- fesim(N=nobs, time=tobs, R=reps, rob=T, crs=crs, tol=1e-10, seed=1565761,
                     b1=1,  bx1a1=0.5,
                     cov_a2_bx1w=p1,
                     bx1w_m=0.5, bx1w_sd=(s1^0.5),
                     a1_m=1, a1_sd=1, a2_m=0.1, a2_sd=(4^0.5),
                     w_m=1, w_sd=(s2^0.5), x1_m=1, x1_sd=(s3^0.5), u_sd=1)
        
        bla <- paste("b",1,s2,s3,sep="_") 
        bl[[bla]][i,j] <- mean(sim$beta_re.df[,1])
        bla <- paste("b",2,s2,s3,sep="_") 
        bl[[bla]][i,j] <- mean(sim$beta_fe.df[,1])
        bla <- paste("b",3,s2,s3,sep="_") 
        bl[[bla]][i,j] <- mean(sim$beta_feis.df[,1])
        bla <- paste("bias",1,s2,s3,sep="_")
        biasl[[bla]][i,j] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
        bla <- paste("bias",2,s2,s3,sep="_")
        biasl[[bla]][i,j] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
        bla <- paste("bias",3,s2,s3,sep="_")
        biasl[[bla]][i,j] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
        bla <- paste("biassd",1,s2,s3,sep="_")
        biassdl[[bla]][i,j] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
        bla <- paste("biassd",2,s2,s3,sep="_")
        biassdl[[bla]][i,j] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
        bla <- paste("biassd",3,s2,s3,sep="_")
        biassdl[[bla]][i,j] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
        bla <- paste("rej",1,s2,s3,sep="_")
        rejl[[bla]][i,j] <- mean(sim$test_p.df[,1]<0.05)
        bla <- paste("rej",2,s2,s3,sep="_")
        rejl[[bla]][i,j] <- mean(sim$test_p.df[,2]<0.05)
        bla <- paste("rej",3,s2,s3,sep="_")
        rejl[[bla]][i,j] <- mean(sim$test_p.df[,3]<0.05)
        bla <- paste("rmse",1,s2,s3,sep="_")
        rmsel[[bla]][i,j] <- rmse(sim$beta_re.df[,1],sim$theta[1])
        bla <- paste("rmse",2,s2,s3,sep="_")
        rmsel[[bla]][i,j] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
        bla <- paste("rmse",3,s2,s3,sep="_")
        rmsel[[bla]][i,j] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
        bla <- paste("bias_pred",2,s2,s3,sep="_")
        bias_predl[[bla]][i,j] <- mean(sim$param_emp.df[,5])
        
        
        save(sim,
             file=paste0("Simulation3_rob_", i, "_", j, "_", k, "_", l, "_smallN.RData", sep=""))
        
        l=l+1
      }
      k=k+1
    }
    j=j+1
  }
  i=i+1
}


### Load data (if using saved simulations)

# i=1
# for (p1 in phi1) {
#   j=1
#   for (s1 in sigma1) {
#     k=1
#     for (s2 in sigma2) {
#       l=1
#       for (s3 in sigma3) {
#         
#         load(paste0("Simulation3_rob_", i, "_", j, "_", k, "_", l, "_smallN.RData", sep=""))
#         
#         bla <- paste("b",1,s2,s3,sep="_") 
#         bl[[bla]][i,j] <- mean(sim$beta_re.df[,1])
#         bla <- paste("b",2,s2,s3,sep="_") 
#         bl[[bla]][i,j] <- mean(sim$beta_fe.df[,1])
#         bla <- paste("b",3,s2,s3,sep="_") 
#         bl[[bla]][i,j] <- mean(sim$beta_feis.df[,1])
#         bla <- paste("bias",1,s2,s3,sep="_")
#         biasl[[bla]][i,j] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
#         bla <- paste("bias",2,s2,s3,sep="_")
#         biasl[[bla]][i,j] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
#         bla <- paste("bias",3,s2,s3,sep="_")
#         biasl[[bla]][i,j] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
#         bla <- paste("biassd",1,s2,s3,sep="_")
#         biassdl[[bla]][i,j] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
#         bla <- paste("biassd",2,s2,s3,sep="_")
#         biassdl[[bla]][i,j] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
#         bla <- paste("biassd",3,s2,s3,sep="_")
#         biassdl[[bla]][i,j] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
#         bla <- paste("rej",1,s2,s3,sep="_")
#         rejl[[bla]][i,j] <- mean(sim$test_p.df[,1]<0.05)
#         bla <- paste("rej",2,s2,s3,sep="_")
#         rejl[[bla]][i,j] <- mean(sim$test_p.df[,2]<0.05)
#         bla <- paste("rej",3,s2,s3,sep="_")
#         rejl[[bla]][i,j] <- mean(sim$test_p.df[,3]<0.05)
#         bla <- paste("rmse",1,s2,s3,sep="_")
#         rmsel[[bla]][i,j] <- rmse(sim$beta_re.df[,1],sim$theta[1])
#         bla <- paste("rmse",2,s2,s3,sep="_")
#         rmsel[[bla]][i,j] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
#         bla <- paste("rmse",3,s2,s3,sep="_")
#         rmsel[[bla]][i,j] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
#         bla <- paste("bias_pred",2,s2,s3,sep="_")
#         bias_predl[[bla]][i,j] <- mean(sim$param_emp.df[,5])
#         
#         l=l+1
#       }
#       k=k+1
#     }
#     j=j+1
#   }
#   i=i+1
# }





#########################
#### Plot results 3) ####
#########################

### Bias: 3D Wireframe Plot

g <- expand.grid(x = phi1, y = sigma1, gr = 1:3)


g$v1 <- as.vector(biasl[["bias_2_0.2_0.2"]])
g$v2 <- as.vector(biasl[["bias_2_0.2_1"]])
g$v3 <- as.vector(biasl[["bias_2_0.2_5"]])
g$v4 <- as.vector(biasl[["bias_2_1_0.2"]])
g$v5 <- as.vector(biasl[["bias_2_1_1"]])
g$v6 <- as.vector(biasl[["bias_2_1_5"]])
g$v7 <- as.vector(biasl[["bias_2_5_0.2"]])
g$v8 <- as.vector(biasl[["bias_2_5_1"]])
g$v9 <- as.vector(biasl[["bias_2_5_5"]])


#attributes(g$v1)
#names(g$v1) <- "bka"

mytitles <- c(expression(paste("G) ", sigma[w]," = 5", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("H) ", sigma[w]," = 5", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("I) ", sigma[w]," = 5", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("D) ", sigma[w]," = 1", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("E) ", sigma[w]," = 1", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("F) ", sigma[w]," = 1", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("A) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("B) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("C) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 5", sep="")))

sty <- list()
sty$axis.line$col <- NA
sty$strip.border$col <- NA
sty$strip.background$col <- NA

limlst <- list(c(0,3), c(0,2), c(0,1), c(0,2), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1))


pl<-wireframe(v7 + v8 + v9 + v4 + v5 + v6 + v1 + v2 + v3 ~ x * y, data = g, outer=T,
              # main="Bias", 
              drape=T, zlab="", 
              xlab = list(label = expression(phi), cex = 1.2), 
              ylab = list(label = expression(sigma[delta]), cex=1.2),
              column.values = list(label = sigma1), 
              row.values = list(label = phi1),  
              strip = strip.custom(factor.levels = mytitles, bg=NA, fg=NA),
              par.strip.text = list(cex=1.5),  
              par.settings = sty,
              scales=list(arrows=F, distance=1.5, 
                          x = list(cex=1.2), y = list(cex=1.2), 
                          z = list(cex=1.2, axs="i")),
              layout = c(3, 3), 
              col.regions=rainbow(75),
              colorkey = list(labels=list(cex=1.5)),
              #screen = list(z = 30, x = -60),
              zlim=c(0,1.5)
)
pl

cairo_pdf(file=paste("../03_Output/Output_", "Simplot3_rob_bias_smallN.pdf", sep=""), width=12.40, height=12.40, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
pl
dev.off()


### Rejection rate: 3D Wireframe Plot

g2 <- expand.grid(x = phi1, y = sigma1, gr = 1:3)


g2$v1 <- as.vector(rejl[["rej_1_0.2_0.2"]])
g2$v2 <- as.vector(rejl[["rej_1_0.2_1"]])
g2$v3 <- as.vector(rejl[["rej_1_0.2_5"]])
g2$v4 <- as.vector(rejl[["rej_1_1_0.2"]])
g2$v5 <- as.vector(rejl[["rej_1_1_1"]])
g2$v6 <- as.vector(rejl[["rej_1_1_5"]])
g2$v7 <- as.vector(rejl[["rej_1_5_0.2"]])
g2$v8 <- as.vector(rejl[["rej_1_5_1"]])
g2$v9 <- as.vector(rejl[["rej_1_5_5"]])
#attributes(g$v1)
#names(g$v1) <- "bka"

mytitles <- c(expression(paste("G) ", sigma[w]," = 5", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("H) ", sigma[w]," = 5", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("I) ", sigma[w]," = 5", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("D) ", sigma[w]," = 1", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("E) ", sigma[w]," = 1", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("F) ", sigma[w]," = 1", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("A) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("B) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("C) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 5", sep="")))

sty <- list()
sty$axis.line$col <- NA
sty$strip.border$col <- NA
sty$strip.background$col <- NA



pl<-wireframe(v7 + v8 + v9 + v4 + v5 + v6 + v1 + v2 + v3 ~ x * y, data = g2, outer=T,
              # main="Bias", 
              drape=T, zlab="", 
              xlab = list(label = expression(phi), cex = 1.2), 
              ylab = list(label = expression(sigma[delta]), cex=1.2),
              column.values = list(label = sigma1), 
              row.values = list(label = phi1),  
              strip = strip.custom(factor.levels = mytitles, bg=NA, fg=NA),
              par.strip.text = list(cex=1.5),  
              par.settings = sty,
              scales=list(arrows=F, distance=1.5, 
                          x = list(cex=1.2), y = list(cex=1.2), z = list(cex=1.2)),
              layout = c(3, 3), 
              col.regions=rainbow(75),
              colorkey = list(labels=list(cex=1.5)),
              #screen = list(z = 30, x = -60),
              zlim=c(0,1)
)
pl

cairo_pdf(file=paste("../03_Output/Output_", "Simplot3_rob_rej_smallN.pdf", sep=""), width=12.40, height=12.40, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
pl
dev.off()




### Rejection rate as function of bias

g3<-Map(cbind, bias=g, rej=g2)
g3<-lapply(g3, FUN=function(x) cbind(data.frame(x), s1=g3$y[,1]))
g3<-g3[4:12]


mytitles <- c(expression(paste("A) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("B) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("C) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("D) ", sigma[w]," = 1", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("E) ", sigma[w]," = 1", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("F) ", sigma[w]," = 1", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("G) ", sigma[w]," = 5", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("H) ", sigma[w]," = 5", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("I) ", sigma[w]," = 5", " ; ", sigma[v]," = 5", sep=""))
)

# Single ggplot

for(i in 1:9){
  pl2 <- ggplot(g3[[i]], aes(x = bias, y = rej)) 
  pl2 <- pl2 + geom_point(aes(x = bias, y = rej, color=as.factor(s1), shape=as.factor(s1)), size=3) 
  pl2 <- pl2 + geom_smooth(aes(y=rej, x=bias), span = 0.5,  se=F, colour="dodgerblue3") 
  pl2 <- pl2 + scale_shape_manual(values=c(3, 15, 16, 17, 18), guide = FALSE)
  pl2 <- pl2 + geom_abline(intercept = 0.05, slope = 0,color="red", linetype=2) 
  pl2 <- pl2 + geom_abline(intercept = 0.95, slope = 0,color="red", linetype=2)
  pl2 <- pl2 + labs(y="", x = "", colour=expression(paste(sigma[delta])), title=mytitles[i]) 
  pl2 <- pl2 + theme_bw() + ylim(0,1.05) # + xlim(0, 0.8)
  pl2 <- pl2 + geom_blank(aes(y = 0)) + geom_blank(aes(y = 1)) # + geom_blank(aes(x = 0.8)) # Fake values for axes
  pl2 <- pl2 + guides(colour = guide_legend(override.aes = list(shape = c(3, 15, 16, 17, 18)) ))
  pl2 <- pl2 + scale_linetype(guide = FALSE)
  pl2 <- pl2 + theme(legend.text = element_text(size = 20),
                     #legend.position="bottom",
                     legend.key=element_blank(),
                     legend.title = element_text(size=20),
                     axis.text.x = element_text(size=16),
                     axis.text.y = element_text(size=16, colour="black"),
                     axis.title.x = element_text(size=20),
                     axis.title.y = element_text(size=20),
                     plot.title = element_text(size=20 ),
                     strip.background =element_blank(),
                     strip.text = element_text(size=20, colour="black"),
                     panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
  
  assign(paste("pl2_", i, sep=""), pl2)
}


# Add axis title to the outer 

figure <- ggarrange(pl2_1, pl2_2, pl2_3, pl2_4, pl2_5, pl2_6, pl2_7, pl2_8, pl2_9,
                    ncol=3, nrow=3, common.legend = TRUE, legend="right",
                    label.x="Bias in FE", label.y="Rejection rate")

figure <- annotate_figure(figure,
                          bottom = text_grob("Bias in FE", size = 20),
                          left = text_grob("Rejection rate", size = 20, rot=90)
)
figure



cairo_pdf(file=paste("../03_Output/Output_", "Simplot4_rob_smallN.pdf", sep=""), width=17.54, height=12.40, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
figure
dev.off()






#################################################################
#### 4) Simulations varying phi, sigma1, sigma2, sigma3: BST ####
#################################################################


# Vary phi1 (Cov alpha2, delta) and sigma1 (Var delta) simultaneously : 
phi1 <- c(seq(0,0.2,0.05),seq(0.3,0.5,0.1),0.65,0.8)
sigma1 <- c(0.2,0.5,1,3,5)

# Vary sigma2 (Var of w) and sigma3 (Var of x)
sigma2 <- c(0.2,1,5)
sigma3 <- c(0.2,1,5)



bl <- list()
biasl <- list()
biassdl <- list()
rejl <- list()
rmsel <- list()
bias_predl <- list()
for (s2 in sigma2) {
  for (s3 in sigma3) {
    for (x in c(1,2,3)) {
      bla <- paste("b",x,s2,s3,sep="_") 
      bl[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
      bla <- paste("bias",x,s2,s3,sep="_") 
      biasl[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
      bla <- paste("biassd",x,s2,s3,sep="_") 
      biassdl[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
      bla <- paste("rej",x,s2,s3,sep="_") 
      rejl[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
      bla <- paste("rmse",x,s2,s3,sep="_") 
      rmsel[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
      bla <- paste("bias_pred",x,s2,s3,sep="_") 
      bias_predl[[bla]] <- matrix(NA, nrow=length(phi1), ncol=length(sigma1))
    }
  }    
}  




### Run simulations

i=1
for (p1 in phi1) {
  j=1
  for (s1 in sigma1) {
    k=1
    for (s2 in sigma2) {
      l=1
      for (s3 in sigma3) {
        
        sim <- fesim_bs(N=nobs, time=tobs, R=reps, bsR=bsreps, crs=crs, 
                        tol=1e-10, seed=1565761,
                        b1=1,  bx1a1=0.5,
                        cov_a2_bx1w=p1,
                        bx1w_m=0.5, bx1w_sd=(s1^0.5),
                        a1_m=1, a1_sd=1, a2_m=0.1, a2_sd=(4^0.5),
                        w_m=1, w_sd=(s2^0.5), x1_m=1, x1_sd=(s3^0.5), u_sd=1)
        
        bla <- paste("b",1,s2,s3,sep="_") 
        bl[[bla]][i,j] <- mean(sim$beta_re.df[,1])
        bla <- paste("b",2,s2,s3,sep="_") 
        bl[[bla]][i,j] <- mean(sim$beta_fe.df[,1])
        bla <- paste("b",3,s2,s3,sep="_") 
        bl[[bla]][i,j] <- mean(sim$beta_feis.df[,1])
        bla <- paste("bias",1,s2,s3,sep="_")
        biasl[[bla]][i,j] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
        bla <- paste("bias",2,s2,s3,sep="_")
        biasl[[bla]][i,j] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
        bla <- paste("bias",3,s2,s3,sep="_")
        biasl[[bla]][i,j] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
        bla <- paste("biassd",1,s2,s3,sep="_")
        biassdl[[bla]][i,j] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
        bla <- paste("biassd",2,s2,s3,sep="_")
        biassdl[[bla]][i,j] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
        bla <- paste("biassd",3,s2,s3,sep="_")
        biassdl[[bla]][i,j] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
        bla <- paste("rej",1,s2,s3,sep="_")
        rejl[[bla]][i,j] <- mean(sim$test_p.df[,1]<0.05)
        bla <- paste("rej",2,s2,s3,sep="_")
        rejl[[bla]][i,j] <- mean(sim$test_p.df[,2]<0.05)
        bla <- paste("rej",3,s2,s3,sep="_")
        rejl[[bla]][i,j] <- mean(sim$test_p.df[,3]<0.05)
        bla <- paste("rmse",1,s2,s3,sep="_")
        rmsel[[bla]][i,j] <- rmse(sim$beta_re.df[,1],sim$theta[1])
        bla <- paste("rmse",2,s2,s3,sep="_")
        rmsel[[bla]][i,j] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
        bla <- paste("rmse",3,s2,s3,sep="_")
        rmsel[[bla]][i,j] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
        bla <- paste("bias_pred",2,s2,s3,sep="_")
        bias_predl[[bla]][i,j] <- mean(sim$param_emp.df[,5])
        
        save(sim,
             file=paste0("Simulation3_bs_", i, "_", j, "_", k, "_", l, "_smallN.RData", sep=""))
        
        l=l+1
      }
      k=k+1
    }
    j=j+1
  }
  i=i+1
}


### Load data (if using saved simulations)

# i=1
# for (p1 in phi1) {
#   j=1
#   for (s1 in sigma1) {
#     k=1
#     for (s2 in sigma2) {
#       l=1
#       for (s3 in sigma3) {
#         
#         load(paste0(wd, "/Simulation3_bs_", i, "_", j, "_", k, "_", l, "_smallN.RData", sep=""))
#         
#         bla <- paste("b",1,s2,s3,sep="_") 
#         bl[[bla]][i,j] <- mean(sim$beta_re.df[,1])
#         bla <- paste("b",2,s2,s3,sep="_") 
#         bl[[bla]][i,j] <- mean(sim$beta_fe.df[,1])
#         bla <- paste("b",3,s2,s3,sep="_") 
#         bl[[bla]][i,j] <- mean(sim$beta_feis.df[,1])
#         bla <- paste("bias",1,s2,s3,sep="_")
#         biasl[[bla]][i,j] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
#         bla <- paste("bias",2,s2,s3,sep="_")
#         biasl[[bla]][i,j] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
#         bla <- paste("bias",3,s2,s3,sep="_")
#         biasl[[bla]][i,j] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
#         bla <- paste("biassd",1,s2,s3,sep="_")
#         biassdl[[bla]][i,j] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
#         bla <- paste("biassd",2,s2,s3,sep="_")
#         biassdl[[bla]][i,j] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
#         bla <- paste("biassd",3,s2,s3,sep="_")
#         biassdl[[bla]][i,j] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
#         bla <- paste("rej",1,s2,s3,sep="_")
#         rejl[[bla]][i,j] <- mean(sim$test_p.df[,1]<0.05)
#         bla <- paste("rej",2,s2,s3,sep="_")
#         rejl[[bla]][i,j] <- mean(sim$test_p.df[,2]<0.05)
#         bla <- paste("rej",3,s2,s3,sep="_")
#         rejl[[bla]][i,j] <- mean(sim$test_p.df[,3]<0.05)
#         bla <- paste("rmse",1,s2,s3,sep="_")
#         rmsel[[bla]][i,j] <- rmse(sim$beta_re.df[,1],sim$theta[1])
#         bla <- paste("rmse",2,s2,s3,sep="_")
#         rmsel[[bla]][i,j] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
#         bla <- paste("rmse",3,s2,s3,sep="_")
#         rmsel[[bla]][i,j] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
#         bla <- paste("bias_pred",2,s2,s3,sep="_")
#         bias_predl[[bla]][i,j] <- mean(sim$param_emp.df[,5])
#         
#         l=l+1
#       }
#       k=k+1
#     }
#     j=j+1
#   }
#   i=i+1
# }





#########################
#### Plot results 4) ####
#########################

### Bias: 3D Wireframe Plot

g <- expand.grid(x = phi1, y = sigma1, gr = 1:3)


g$v1 <- as.vector(biasl[["bias_2_0.2_0.2"]])
g$v2 <- as.vector(biasl[["bias_2_0.2_1"]])
g$v3 <- as.vector(biasl[["bias_2_0.2_5"]])
g$v4 <- as.vector(biasl[["bias_2_1_0.2"]])
g$v5 <- as.vector(biasl[["bias_2_1_1"]])
g$v6 <- as.vector(biasl[["bias_2_1_5"]])
g$v7 <- as.vector(biasl[["bias_2_5_0.2"]])
g$v8 <- as.vector(biasl[["bias_2_5_1"]])
g$v9 <- as.vector(biasl[["bias_2_5_5"]])


#attributes(g$v1)
#names(g$v1) <- "bka"

mytitles <- c(expression(paste("G) ", sigma[w]," = 5", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("H) ", sigma[w]," = 5", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("I) ", sigma[w]," = 5", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("D) ", sigma[w]," = 1", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("E) ", sigma[w]," = 1", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("F) ", sigma[w]," = 1", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("A) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("B) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("C) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 5", sep="")))

sty <- list()
sty$axis.line$col <- NA
sty$strip.border$col <- NA
sty$strip.background$col <- NA

limlst <- list(c(0,3), c(0,2), c(0,1), c(0,2), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1))



pl<-wireframe(v7 + v8 + v9 + v4 + v5 + v6 + v1 + v2 + v3 ~ x * y, data = g, outer=T,
              # main="Bias", 
              drape=T, zlab="", 
              xlab = list(label = expression(phi), cex = 1.2), 
              ylab = list(label = expression(sigma[delta]), cex=1.2),
              column.values = list(label = sigma1), 
              row.values = list(label = phi1),  
              strip = strip.custom(factor.levels = mytitles, bg=NA, fg=NA),
              par.strip.text = list(cex=1.5),  
              par.settings = sty,
              scales=list(arrows=F, distance=1.5, 
                          x = list(cex=1.2), y = list(cex=1.2), 
                          z = list(cex=1.2, axs="i")),
              layout = c(3, 3), 
              col.regions=rainbow(75),
              colorkey = list(labels=list(cex=1.5)),
              #screen = list(z = 30, x = -60),
              zlim=c(0,1.5)
)
pl

cairo_pdf(file=paste("../03_Output/Output_", "Simplot3_bs_bias_smallN.pdf", sep=""), width=12.40, height=12.40, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
pl
dev.off()




bias_predl$bias_pred_2_5_0.1
biasl$bias_2_5_0.1



### Rejection rate: 3D Wireframe Plot

g2 <- expand.grid(x = phi1, y = sigma1, gr = 1:3)


g2$v1 <- as.vector(rejl[["rej_1_0.2_0.2"]])
g2$v2 <- as.vector(rejl[["rej_1_0.2_1"]])
g2$v3 <- as.vector(rejl[["rej_1_0.2_5"]])
g2$v4 <- as.vector(rejl[["rej_1_1_0.2"]])
g2$v5 <- as.vector(rejl[["rej_1_1_1"]])
g2$v6 <- as.vector(rejl[["rej_1_1_5"]])
g2$v7 <- as.vector(rejl[["rej_1_5_0.2"]])
g2$v8 <- as.vector(rejl[["rej_1_5_1"]])
g2$v9 <- as.vector(rejl[["rej_1_5_5"]])
#attributes(g$v1)
#names(g$v1) <- "bka"

mytitles <- c(expression(paste("G) ", sigma[w]," = 5", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("H) ", sigma[w]," = 5", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("I) ", sigma[w]," = 5", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("D) ", sigma[w]," = 1", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("E) ", sigma[w]," = 1", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("F) ", sigma[w]," = 1", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("A) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("B) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("C) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 5", sep="")))

sty <- list()
sty$axis.line$col <- NA
sty$strip.border$col <- NA
sty$strip.background$col <- NA



pl<-wireframe(v7 + v8 + v9 + v4 + v5 + v6 + v1 + v2 + v3 ~ x * y, data = g2, outer=T,
              # main="Bias", 
              drape=T, zlab="", 
              xlab = list(label = expression(phi), cex = 1.2), 
              ylab = list(label = expression(sigma[delta]), cex=1.2),
              column.values = list(label = sigma1), 
              row.values = list(label = phi1),  
              strip = strip.custom(factor.levels = mytitles, bg=NA, fg=NA),
              par.strip.text = list(cex=1.5),  
              par.settings = sty,
              scales=list(arrows=F, distance=1.5, 
                          x = list(cex=1.2), y = list(cex=1.2), z = list(cex=1.2)),
              layout = c(3, 3), 
              col.regions=rainbow(75),
              colorkey = list(labels=list(cex=1.5)),
              #screen = list(z = 30, x = -60),
              zlim=c(0,1)
)
pl

cairo_pdf(file=paste("../03_Output/Output_", "Simplot3_bs_rej_smallN.pdf", sep=""), width=12.40, height=12.40, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
pl
dev.off()





### Rejection rate as function of bias

g3<-Map(cbind, bias=g, rej=g2)
g3<-lapply(g3, FUN=function(x) cbind(data.frame(x), s1=g3$y[,1]))
g3<-g3[4:12]


mytitles <- c(expression(paste("A) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("B) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("C) ", sigma[w]," = 0.2", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("D) ", sigma[w]," = 1", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("E) ", sigma[w]," = 1", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("F) ", sigma[w]," = 1", " ; ", sigma[v]," = 5", sep="")),
              expression(paste("G) ", sigma[w]," = 5", " ; ", sigma[v]," = 0.2", sep="")),
              expression(paste("H) ", sigma[w]," = 5", " ; ", sigma[v]," = 1", sep="")),
              expression(paste("I) ", sigma[w]," = 5", " ; ", sigma[v]," = 5", sep=""))
)

# Single ggplot

for(i in 1:9){
  pl2 <- ggplot(g3[[i]], aes(x = bias, y = rej)) 
  pl2 <- pl2 + geom_point(aes(x = bias, y = rej, color=as.factor(s1), shape=as.factor(s1)), size=3) 
  pl2 <- pl2 + geom_smooth(aes(y=rej, x=bias), span = 0.5,  se=F, colour="dodgerblue3") 
  pl2 <- pl2 + scale_shape_manual(values=c(3, 15, 16, 17, 18), guide = FALSE)
  pl2 <- pl2 + geom_abline(intercept = 0.05, slope = 0,color="red", linetype=2) 
  pl2 <- pl2 + geom_abline(intercept = 0.95, slope = 0,color="red", linetype=2)
  pl2 <- pl2 + labs(y="", x = "", colour=expression(paste(sigma[delta])), title=mytitles[i]) 
  pl2 <- pl2 + theme_bw() + ylim(0,1.05) # + xlim(0, 0.8)
  pl2 <- pl2 + geom_blank(aes(y = 0)) + geom_blank(aes(y = 1)) # + geom_blank(aes(x = 0.8)) # Fake values for axes
  pl2 <- pl2 + guides(colour = guide_legend(override.aes = list(shape = c(3, 15, 16, 17, 18)) ))
  pl2 <- pl2 + scale_linetype(guide = FALSE)
  pl2 <- pl2 + theme(legend.text = element_text(size = 20),
                     #legend.position="bottom",
                     legend.key=element_blank(),
                     legend.title = element_text(size=20),
                     axis.text.x = element_text(size=16),
                     axis.text.y = element_text(size=16, colour="black"),
                     axis.title.x = element_text(size=20),
                     axis.title.y = element_text(size=20),
                     plot.title = element_text(size=20 ),
                     strip.background =element_blank(),
                     strip.text = element_text(size=20, colour="black"),
                     panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
  
  assign(paste("pl2_", i, sep=""), pl2)
}


# Add axis title to the outer 

figure <- ggarrange(pl2_1, pl2_2, pl2_3, pl2_4, pl2_5, pl2_6, pl2_7, pl2_8, pl2_9,
                    ncol=3, nrow=3, common.legend = TRUE, legend="right",
                    label.x="Bias in FE", label.y="Rejection rate")

figure <- annotate_figure(figure,
                          bottom = text_grob("Bias in FE", size = 20),
                          left = text_grob("Rejection rate", size = 20, rot=90)
)
figure



cairo_pdf(file=paste("../03_Output/Output_", "Simplot4_bs_smallN.pdf", sep=""), width=17.54, height=12.40, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
figure
dev.off()














####----------------------------------------------------------------------####
####                   Time lags, Sample selection                        ####
####----------------------------------------------------------------------####


############################
#### Specify parameters ####
############################

### Number of obs per sim
nobs <- 300

### Time periods per sim
tobs <- 10

### Replications per sim
reps <- 1000

### Replications per bootstrap run
bsreps <- 100

### Number of CPU cores used (for HPC verion of program functions)
# If set to NA in combination with hpc program verion, all cores are used (only recommended for hpc)
crs <- NA





#############################################################
#### 5) Use time lags: Simulations varying phi, lag: ART ####
#############################################################


### Varying phi1:
# covariance a2 and phi1 (parameter: cov_a2_bx1w)
phi1 <- c(-0.8,seq(-0.5,-0.3,0.1),seq(-0.25,0.25,0.05),seq(0.3,0.5,0.1),0.8)

### Varying proportion of lagged effects on bxt effect
lags <- c(0.1, 0.3, 0.5, 1) 

### results matrix
df0_4 <- data.frame(phi1=phi1)
df1_4 <- df0_4
df2_4 <- df0_4
df3_4 <- df0_4


start_time <- Sys.time()
### Simulations
for (l in lags) {
  i=1
  for (p1 in phi1) {
    sim <- fesim(N=nobs, time=tobs, R=reps, rob=T, crs=crs, tol=1e-10, seed=1565761, lag = TRUE,
                 b1=1, bx1a1=1, lagprop=l, 
                 cov_a2_bx1w=p1,
                 bx1w_m=1, bx1w_sd=(1^0.5),
                 a1_m=1, a1_sd=(4^0.5), a2_m=0.1, a2_sd=(4^0.5),
                 w_m=1, w_sd=1, x1_m=1, x1_sd=(1^0.5), u_sd=1)
    if (l==lags[1]) {
      #relative bias
      df0_4$b1[i] <- mean(sim$beta_re.df[,1])
      df0_4$b2[i] <- mean(sim$beta_fe.df[,1])
      df0_4$b3[i] <- mean(sim$beta_feis.df[,1])
      df0_4$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df0_4$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df0_4$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df0_4$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0_4$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0_4$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0_4$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df0_4$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df0_4$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df0_4$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df0_4$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df0_4$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df0_4$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      df0_4$theta_re[i] <- mean(sim$theta_re.df)
      df0_4$Va_re[i] <- mean(sim$Va_re.df)
      df0_4$Ve_re[i] <- mean(sim$Ve_re.df)
      df0_4$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    if (l==lags[2]) {
      #relative bias
      df1_4$b1[i] <- mean(sim$beta_re.df[,1])
      df1_4$b2[i] <- mean(sim$beta_fe.df[,1])
      df1_4$b3[i] <- mean(sim$beta_feis.df[,1])
      df1_4$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df1_4$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df1_4$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df1_4$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1_4$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1_4$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1_4$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df1_4$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df1_4$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df1_4$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df1_4$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df1_4$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df1_4$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      df1_4$theta_re[i] <- mean(sim$theta_re.df)
      df1_4$Va_re[i] <- mean(sim$Va_re.df)
      df1_4$Ve_re[i] <- mean(sim$Ve_re.df)
      df1_4$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    if (l==lags[3]) {
      #relative bias
      df2_4$b1[i] <- mean(sim$beta_re.df[,1])
      df2_4$b2[i] <- mean(sim$beta_fe.df[,1])
      df2_4$b3[i] <- mean(sim$beta_feis.df[,1])
      df2_4$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df2_4$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df2_4$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df2_4$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df2_4$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df2_4$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df2_4$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df2_4$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df2_4$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df2_4$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df2_4$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df2_4$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df2_4$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      df2_4$theta_re[i] <- mean(sim$theta_re.df)
      df2_4$Va_re[i] <- mean(sim$Va_re.df)
      df2_4$Ve_re[i] <- mean(sim$Ve_re.df)
      df2_4$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    if (l==lags[4]) {
      #relative bias
      df3_4$b1[i] <- mean(sim$beta_re.df[,1])
      df3_4$b2[i] <- mean(sim$beta_fe.df[,1])
      df3_4$b3[i] <- mean(sim$beta_feis.df[,1])
      df3_4$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df3_4$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df3_4$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df3_4$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df3_4$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df3_4$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df3_4$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df3_4$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df3_4$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df3_4$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df3_4$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df3_4$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df3_4$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      df3_4$theta_re[i] <- mean(sim$theta_re.df)
      df3_4$Va_re[i] <- mean(sim$Va_re.df)
      df3_4$Ve_re[i] <- mean(sim$Ve_re.df)
      df3_4$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    
    save(sim,
         file=paste0("Simulation4_", which(lags == l), "_", i, ".RData", sep=""))
    
    i=i+1
  }
}
end_time <- Sys.time()
end_time - start_time

save(df0_4, df1_4, df2_4, df3_4, file="ART4.RData")




################################################
### Plot results 5)                         ####
################################################






### 1) results for prolag 1
res <- reshape(df0_4, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 

library(ggplot2)
res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot4_1 <- ggplot(res, aes(phi1, bias1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste(beta[xn-1], " = ", beta[xn-2], " = ", v, beta[xn]), list(v=lags[1]))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.45,1.2), breaks=c(seq(-0.4,-0.2,0.2), 0, seq(0.2,1.2,0.2))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot4_1


### 2) results for prolag 2
res <- reshape(df1_4, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 

library(ggplot2)
res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot4_2 <- ggplot(res, aes(phi1, bias1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste(beta[xn-1], " = ", beta[xn-2], " = ", v, beta[xn]), list(v=lags[2]))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.45,1.2), breaks=c(seq(-0.4,-0.2,0.2), 0, seq(0.2,1.2,0.2))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot4_2


### 3) results for prolag 3
res <- reshape(df2_4, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 

library(ggplot2)
res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot4_3 <- ggplot(res, aes(phi1, bias1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste(beta[xn-1], " = ", beta[xn-2], " = ", v, beta[xn]), list(v=lags[3]))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.45,1.2), breaks=c(seq(-0.4,-0.2,0.2), 0, seq(0.2,1.2,0.2))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot4_3


### 4) results for prolag 4
res <- reshape(df3_4, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 

library(ggplot2)
res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot4_4 <- ggplot(res, aes(phi1, bias1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste(beta[xn-1], " = ", beta[xn-2], " = ", v, beta[xn]), list(v=lags[4]))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.45,1.2), breaks=c(seq(-0.4,-0.2,0.2), 0, seq(0.2,1.2,0.2))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.9,0.3), legend.justification=c(0.9,0.3),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot4_4

### Combine Figures
library(ggpubr)
Simplot4 <- ggarrange(simplot4_1,simplot4_2,simplot4_3,simplot4_4,
                      ncol=2, nrow=2, common.legend = F, 
                      labels=c("A) ", "B) ", "C) ", "D) "), font.label = list(size = 20)) 

cairo_pdf(file=paste("../03_Output/Output_", "Simplot5_lag.pdf", sep=""), width=17.54, height=10, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
Simplot4
dev.off()


##### Art Rejection rate


### 1) Art rejection for prolag 1
res <- reshape(df0_4, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot5_1 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste(beta[xn-1], " = ", beta[xn-2], " = ", v, beta[xn]), list(v=lags[1]))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot5_1


### 2) Art rejection for prolag 2
res <- reshape(df1_4, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot5_2 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste(beta[xn-1], " = ", beta[xn-2], " = ", v, beta[xn]), list(v=lags[2]))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot5_2


### 3) Art rejection for prolag 3
res <- reshape(df2_4, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot5_3 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste(beta[xn-1], " = ", beta[xn-2], " = ", v, beta[xn]), list(v=lags[3]))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot5_3


### 4) Art rejection for prolag 4
res <- reshape(df3_4, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot5_4 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste(beta[xn-1], " = ", beta[xn-2], " = ", v, beta[xn]), list(v=lags[4]))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot5_4



### Combine Figures
library(ggpubr)
Simplot5 <- ggarrange(simplot5_1,simplot5_2,simplot5_3,simplot5_4,
                      ncol=2, nrow=2, common.legend = F, 
                      labels=c("A) ", 
                               "B) ", 
                               "C) ", 
                               "D) "
                      ), font.label = list(size = 20)) 

cairo_pdf(file=paste("../03_Output/Output_", "Simplot5_lag_rejection.pdf", sep=""), width=17.54, height=10, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
Simplot5
dev.off()















################################################################################
#### 6) Panel selection / drop out: Simulations varying phi, selection: ART ####
################################################################################


### Varying phi1:
# covariance a2 and phi1 (parameter: cov_a2_bx1w)
phi1 <- c(-0.8,seq(-0.5,-0.3,0.1),seq(-0.25,0.25,0.05),seq(0.3,0.5,0.1),0.8)

### Varying panel selection effects
sel <- list(list(unbal = 1, selectivity = 0.6),
            list(unbal = 1, selectivity = 0.8),
            list(unbal = 2, selectivity = 0.6),
            list(unbal = 2, selectivity = 0.8)) 

### results matrix
df0_5 <- data.frame(phi1=phi1)
df1_5 <- df0_5
df2_5 <- df0_5
df3_5 <- df0_5


start_time <- Sys.time()
### Simulations
for (l in 1:length(sel)) {
  i=1
  for (p1 in phi1) {
    sim <- fesim(N=nobs, time=tobs, R=reps, rob=T, crs=crs, tol=1e-10, seed=1565761, lag = FALSE,
                 b1=1, bx1a1=1, unbal = sel[[l]]$unbal, selectivity =  sel[[l]]$selectivity,
                 cov_a2_bx1w=p1,
                 bx1w_m=1, bx1w_sd=(1^0.5),
                 a1_m=1, a1_sd=(4^0.5), a2_m=0.1, a2_sd=(4^0.5),
                 w_m=1, w_sd=1, x1_m=1, x1_sd=(1^0.5), u_sd=1)
    if (l==1) {
      #relative bias
      df0_5$b1[i] <- mean(sim$beta_re.df[,1])
      df0_5$b2[i] <- mean(sim$beta_fe.df[,1])
      df0_5$b3[i] <- mean(sim$beta_feis.df[,1])
      df0_5$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df0_5$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df0_5$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df0_5$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0_5$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0_5$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0_5$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df0_5$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df0_5$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df0_5$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df0_5$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df0_5$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df0_5$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      df0_5$theta_re[i] <- mean(sim$theta_re.df)
      df0_5$Va_re[i] <- mean(sim$Va_re.df)
      df0_5$Ve_re[i] <- mean(sim$Ve_re.df)
      df0_5$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    if (l==2) {
      #relative bias
      df1_5$b1[i] <- mean(sim$beta_re.df[,1])
      df1_5$b2[i] <- mean(sim$beta_fe.df[,1])
      df1_5$b3[i] <- mean(sim$beta_feis.df[,1])
      df1_5$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df1_5$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df1_5$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df1_5$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1_5$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1_5$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1_5$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df1_5$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df1_5$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df1_5$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df1_5$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df1_5$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df1_5$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      df1_5$theta_re[i] <- mean(sim$theta_re.df)
      df1_5$Va_re[i] <- mean(sim$Va_re.df)
      df1_5$Ve_re[i] <- mean(sim$Ve_re.df)
      df1_5$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    if (l==3) {
      #relative bias
      df2_5$b1[i] <- mean(sim$beta_re.df[,1])
      df2_5$b2[i] <- mean(sim$beta_fe.df[,1])
      df2_5$b3[i] <- mean(sim$beta_feis.df[,1])
      df2_5$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df2_5$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df2_5$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df2_5$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df2_5$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df2_5$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df2_5$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df2_5$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df2_5$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df2_5$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df2_5$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df2_5$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df2_5$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      df2_5$theta_re[i] <- mean(sim$theta_re.df)
      df2_5$Va_re[i] <- mean(sim$Va_re.df)
      df2_5$Ve_re[i] <- mean(sim$Ve_re.df)
      df2_5$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    if (l==4) {
      #relative bias
      df3_5$b1[i] <- mean(sim$beta_re.df[,1])
      df3_5$b2[i] <- mean(sim$beta_fe.df[,1])
      df3_5$b3[i] <- mean(sim$beta_feis.df[,1])
      df3_5$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df3_5$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df3_5$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df3_5$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df3_5$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df3_5$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df3_5$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df3_5$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df3_5$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      df3_5$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      df3_5$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      df3_5$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      df3_5$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      df3_5$theta_re[i] <- mean(sim$theta_re.df)
      df3_5$Va_re[i] <- mean(sim$Va_re.df)
      df3_5$Ve_re[i] <- mean(sim$Ve_re.df)
      df3_5$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    
    save(sim,
         file=paste0("Simulation5_", l, "_", i, ".RData", sep=""))
    
    i=i+1
  }
}
end_time <- Sys.time()
end_time - start_time

save(df0_5, df1_5, df2_5, df3_5, file="ART5.RData")







################################################
### Plot results 6)                         ####
################################################


# load("ART5.RData")




### 1) results for selection 1
res <- reshape(df0_5, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 

library(ggplot2)
res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot6_1 <- ggplot(res, aes(phi1, bias1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste("Selection on response: ", "p=", v ), list(v=sel[[1]]$selectivity))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.74,0.52), breaks=c(seq(-0.7,-0.1,0.1), 0, seq(0.1,0.5,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot6_1


### 2) results for selection 2
res <- reshape(df1_5, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 

library(ggplot2)
res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot6_2 <- ggplot(res, aes(phi1, bias1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste("Selection on response: ", "p=", v ), list(v=sel[[2]]$selectivity))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.74,0.52), breaks=c(seq(-0.7,-0.1,0.1), 0, seq(0.1,0.5,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot6_2


### 3) results for slection 3
res <- reshape(df2_5, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 

library(ggplot2)
res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot6_3 <- ggplot(res, aes(phi1, bias1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste("Selection on heterogeneous treatment: ", "q=", v ), list(v=sel[[3]]$selectivity))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.5,0.8), breaks=c(seq(-0.5,-0.1,0.1), 0, seq(0.1,0.8,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot6_3


### 4) results for selection 4
res <- reshape(df3_5, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 

library(ggplot2)
res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot6_5 <- ggplot(res, aes(phi1, bias1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste("Selection on heterogeneous treatment: ", "q=", v ), list(v=sel[[4]]$selectivity))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.5,0.8), breaks=c(seq(-0.5,-0.1,0.1), 0, seq(0.1,0.8,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot6_5

### Combine Figures
library(ggpubr)
Simplot6 <- ggarrange(simplot6_1,simplot6_2,simplot6_3,simplot6_5,
                      ncol=2, nrow=2, common.legend = F, 
                      labels=c("A) ", "B) ", "C) ", "D) "), font.label = list(size = 20)) 

cairo_pdf(file=paste("../03_Output/Output_", "Simplot6_sel.pdf", sep=""), width=17.54, height=10, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
Simplot6
dev.off()


##### Art Rejection rate


### 1) Art rejection for selection 1
res <- reshape(df0_5, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot7_1 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste("Selection on response: ", "p=", v ), list(v=sel[[2]]$selectivity))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot7_1


### 2) Art rejection for selection 2
res <- reshape(df1_5, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot7_2 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste("Selection on response: ", "p=", v ), list(v=sel[[2]]$selectivity))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot7_2


### 3) Art rejection for selection 3
res <- reshape(df2_5, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot7_3 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste("Selection on heterogeneous treatment: ", "q=", v ), list(v=sel[[3]]$selectivity))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot7_3


### 4) Art rejection for selection 4
res <- reshape(df3_5, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16))
res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot7_5 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste("Selection on heterogeneous treatment: ", "q=", v ), list(v=sel[[4]]$selectivity))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot7_5



### Combine Figures
library(ggpubr)
Simplot7 <- ggarrange(simplot7_1,simplot7_2,simplot7_3,simplot7_5,
                      ncol=2, nrow=2, common.legend = F, 
                      labels=c("A) ", 
                               "B) ", 
                               "C) ", 
                               "D) "
                      ), font.label = list(size = 20)) 

cairo_pdf(file=paste("../03_Output/Output_", "Simplot6_sel_rejection.pdf", sep=""), width=17.54, height=10, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
Simplot7
dev.off()











###################################################################################
#### 7) Use time lags in sim and estimation: Simulations varying phi, lag: ART ####
###################################################################################


### Varying phi1:
# covariance a2 and phi1 (parameter: cov_a2_bx1w)
phi1 <- c(-0.8,seq(-0.5,-0.3,0.1),seq(-0.25,0.25,0.05),seq(0.3,0.5,0.1),0.8)

### Varying proportion of lagged effects on bxt effect
lags <- c(0.1, 0.3, 0.5, 1) 

### results matrix
df0_6 <- data.frame(phi1=phi1)
df1_6 <- df0_6
df2_6 <- df0_6
df3_6 <- df0_6


start_time <- Sys.time()
### Simulations
for (l in lags) {
  i=1
  for (p1 in phi1) {
    sim <- fesim(N=nobs, time=tobs, R=reps, rob=T, crs=crs, tol=1e-10, seed=1565761, lag = TRUE,
                 b1=1, bx1a1=1, lagprop=l, lagest = TRUE,
                 cov_a2_bx1w=p1,
                 bx1w_m=1, bx1w_sd=(1^0.5),
                 a1_m=1, a1_sd=(4^0.5), a2_m=0.1, a2_sd=(4^0.5),
                 w_m=1, w_sd=1, x1_m=1, x1_sd=(1^0.5), u_sd=1)
    if (l==lags[1]) {
      #relative bias
      df0_6$b1[i] <- mean(sim$beta_re.df[,1])
      df0_6$b2[i] <- mean(sim$beta_fe.df[,1])
      df0_6$b3[i] <- mean(sim$beta_feis.df[,1])
      df0_6$b1l1[i] <- mean(sim$beta_re.df[,2])
      df0_6$b2l1[i] <- mean(sim$beta_fe.df[,2])
      df0_6$b3l1[i] <- mean(sim$beta_feis.df[,2])
      df0_6$b1l2[i] <- mean(sim$beta_re.df[,3])
      df0_6$b2l2[i] <- mean(sim$beta_fe.df[,3])
      df0_6$b3l2[i] <- mean(sim$beta_feis.df[,3])
      
      df0_6$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df0_6$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df0_6$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df0_6$bias1l1[i] <- mean((sim$beta_re.df[,2]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df0_6$bias2l1[i] <- mean((sim$beta_fe.df[,2]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df0_6$bias3l1[i] <- mean((sim$beta_feis.df[,2]- sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df0_6$bias1l2[i] <- mean((sim$beta_re.df[,3]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df0_6$bias2l2[i] <- mean((sim$beta_fe.df[,3]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df0_6$bias3l2[i] <- mean((sim$beta_feis.df[,3]- sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      
      df0_6$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0_6$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0_6$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df0_6$biassd1l1[i] <- var((sim$beta_re.df[,2]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df0_6$biassd2l1[i] <- var((sim$beta_fe.df[,2]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df0_6$biassd3l1[i] <- var((sim$beta_feis.df[,2]-  sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df0_6$biassd1l2[i] <- var((sim$beta_re.df[,3]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df0_6$biassd2l2[i] <- var((sim$beta_fe.df[,3]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df0_6$biassd3l2[i] <- var((sim$beta_feis.df[,3]-  sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      
      df0_6$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df0_6$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df0_6$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
      # df0_6$rmse1[i] <- rmse(sim$beta_re.df[,1],sim$theta[1])
      # df0_6$rmse2[i] <- rmse(sim$beta_fe.df[,1],sim$theta[1])
      # df0_6$rmse3[i] <- rmse(sim$beta_feis.df[,1],sim$theta[1])
      # df0_6$bias2_pred[i] <- mean(sim$param_emp.df[,5]/sim$theta[1])
      # df0_6$theta_re[i] <- mean(sim$theta_re.df)
      # df0_6$Va_re[i] <- mean(sim$Va_re.df)
      # df0_6$Ve_re[i] <- mean(sim$Ve_re.df)
      # df0_6$theta_re2[i] <- mean(1- (sim$Ve_re.df /(10*sim$Va_re.df+sim$Ve_re.df))^.5)
    }
    if (l==lags[2]) {
      #relative bias
      df1_6$b1[i] <- mean(sim$beta_re.df[,1])
      df1_6$b2[i] <- mean(sim$beta_fe.df[,1])
      df1_6$b3[i] <- mean(sim$beta_feis.df[,1])
      df1_6$b1l1[i] <- mean(sim$beta_re.df[,2])
      df1_6$b2l1[i] <- mean(sim$beta_fe.df[,2])
      df1_6$b3l1[i] <- mean(sim$beta_feis.df[,2])
      df1_6$b1l2[i] <- mean(sim$beta_re.df[,3])
      df1_6$b2l2[i] <- mean(sim$beta_fe.df[,3])
      df1_6$b3l2[i] <- mean(sim$beta_feis.df[,3])
      
      df1_6$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df1_6$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df1_6$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df1_6$bias1l1[i] <- mean((sim$beta_re.df[,2]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df1_6$bias2l1[i] <- mean((sim$beta_fe.df[,2]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df1_6$bias3l1[i] <- mean((sim$beta_feis.df[,2]- sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df1_6$bias1l2[i] <- mean((sim$beta_re.df[,3]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df1_6$bias2l2[i] <- mean((sim$beta_fe.df[,3]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df1_6$bias3l2[i] <- mean((sim$beta_feis.df[,3]- sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      
      df1_6$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1_6$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1_6$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df1_6$biassd1l1[i] <- var((sim$beta_re.df[,2]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df1_6$biassd2l1[i] <- var((sim$beta_fe.df[,2]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df1_6$biassd3l1[i] <- var((sim$beta_feis.df[,2]-  sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df1_6$biassd1l2[i] <- var((sim$beta_re.df[,3]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df1_6$biassd2l2[i] <- var((sim$beta_fe.df[,3]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df1_6$biassd3l2[i] <- var((sim$beta_feis.df[,3]-  sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      
      df1_6$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df1_6$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df1_6$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
    }
    if (l==lags[3]) {
      #relative bias
      df2_6$b1[i] <- mean(sim$beta_re.df[,1])
      df2_6$b2[i] <- mean(sim$beta_fe.df[,1])
      df2_6$b3[i] <- mean(sim$beta_feis.df[,1])
      df2_6$b1l1[i] <- mean(sim$beta_re.df[,2])
      df2_6$b2l1[i] <- mean(sim$beta_fe.df[,2])
      df2_6$b3l1[i] <- mean(sim$beta_feis.df[,2])
      df2_6$b1l2[i] <- mean(sim$beta_re.df[,3])
      df2_6$b2l2[i] <- mean(sim$beta_fe.df[,3])
      df2_6$b3l2[i] <- mean(sim$beta_feis.df[,3])
      
      df2_6$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df2_6$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df2_6$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df2_6$bias1l1[i] <- mean((sim$beta_re.df[,2]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df2_6$bias2l1[i] <- mean((sim$beta_fe.df[,2]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df2_6$bias3l1[i] <- mean((sim$beta_feis.df[,2]- sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df2_6$bias1l2[i] <- mean((sim$beta_re.df[,3]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df2_6$bias2l2[i] <- mean((sim$beta_fe.df[,3]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df2_6$bias3l2[i] <- mean((sim$beta_feis.df[,3]- sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      
      df2_6$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df2_6$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df2_6$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df2_6$biassd1l1[i] <- var((sim$beta_re.df[,2]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df2_6$biassd2l1[i] <- var((sim$beta_fe.df[,2]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df2_6$biassd3l1[i] <- var((sim$beta_feis.df[,2]-  sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df2_6$biassd1l2[i] <- var((sim$beta_re.df[,3]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df2_6$biassd2l2[i] <- var((sim$beta_fe.df[,3]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df2_6$biassd3l2[i] <- var((sim$beta_feis.df[,3]-  sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      
      df2_6$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df2_6$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df2_6$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
    }
    if (l==lags[4]) {
      #relative bias
      df3_6$b1[i] <- mean(sim$beta_re.df[,1])
      df3_6$b2[i] <- mean(sim$beta_fe.df[,1])
      df3_6$b3[i] <- mean(sim$beta_feis.df[,1])
      df3_6$b1l1[i] <- mean(sim$beta_re.df[,2])
      df3_6$b2l1[i] <- mean(sim$beta_fe.df[,2])
      df3_6$b3l1[i] <- mean(sim$beta_feis.df[,2])
      df3_6$b1l2[i] <- mean(sim$beta_re.df[,3])
      df3_6$b2l2[i] <- mean(sim$beta_fe.df[,3])
      df3_6$b3l2[i] <- mean(sim$beta_feis.df[,3])
      
      df3_6$bias1[i] <- mean((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])
      df3_6$bias2[i] <- mean((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])
      df3_6$bias3[i] <- mean((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])
      df3_6$bias1l1[i] <- mean((sim$beta_re.df[,2]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df3_6$bias2l1[i] <- mean((sim$beta_fe.df[,2]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df3_6$bias3l1[i] <- mean((sim$beta_feis.df[,2]- sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df3_6$bias1l2[i] <- mean((sim$beta_re.df[,3]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df3_6$bias2l2[i] <- mean((sim$beta_fe.df[,3]-   sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      df3_6$bias3l2[i] <- mean((sim$beta_feis.df[,3]- sim$lagprop * sim$theta[1]) / (sim$lagprop* sim$theta[1]))
      
      df3_6$biassd1[i] <- var((sim$beta_re.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df3_6$biassd2[i] <- var((sim$beta_fe.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df3_6$biassd3[i] <- var((sim$beta_feis.df[,1]-sim$theta[1]) / sim$theta[1])^0.5
      df3_6$biassd1l1[i] <- var((sim$beta_re.df[,2]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df3_6$biassd2l1[i] <- var((sim$beta_fe.df[,2]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df3_6$biassd3l1[i] <- var((sim$beta_feis.df[,2]-  sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df3_6$biassd1l2[i] <- var((sim$beta_re.df[,3]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df3_6$biassd2l2[i] <- var((sim$beta_fe.df[,3]-    sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      df3_6$biassd3l2[i] <- var((sim$beta_feis.df[,3]-  sim$lagprop * sim$theta[1]) / sim$lagprop * sim$theta[1])^0.5
      
      df3_6$rej1[i] <- mean(sim$test_p.df[,1]<0.05)
      df3_6$rej2[i] <- mean(sim$test_p.df[,2]<0.05)
      df3_6$rej3[i] <- mean(sim$test_p.df[,3]<0.05)
    }
    
    save(sim,
         file=paste0("Simulation6_", which(lags == l), "_", i, ".RData", sep=""))
    
    i=i+1
  }
}
end_time <- Sys.time()
end_time - start_time

save(df0_6, df1_6, df2_6, df3_6, file="ART6.RData")





################################################
### Plot results 7)                         ####
################################################


# load("ART6.RData")


#### zeta = 0.5

### 1) results for xt
res <- reshape(df2_6, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16,17:19,20:22,23:25,26:28,29:31))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 

res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot8_1 <- ggplot(res, aes(phi1, bias1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias in % of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste("Bias in ", x[n], "; ", zeta, "=", v), list(v=lags[3]))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.45,0.5), breaks=c(seq(-0.4,-0.1,0.1), 0, seq(0.1,0.5,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot8_1


### 2) results for xt-1
res <- reshape(df2_6, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16,17:19,20:22,23:25,26:28,29:31))
res$up <- res$bias1l1 + 2*res$biassd1l1 
res$lo <- res$bias1l1 - 2*res$biassd1l1 

res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot8_2 <- ggplot(res, aes(phi1, bias1l1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias in % of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste("Bias in ", x[n-1], "; ", zeta, "=", v), list(v=lags[3]))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.45,0.5), breaks=c(seq(-0.4,-0.1,0.1), 0, seq(0.1,0.5,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot8_2


### 3) results for xt-2
res <- reshape(df2_6, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16,17:19,20:22,23:25,26:28,29:31))
res$up <- res$bias1l2 + 2*res$biassd1l2 
res$lo <- res$bias1l2 - 2*res$biassd1l2 

res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot8_3 <- ggplot(res, aes(phi1, bias1l2, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias in % of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste("Bias in ", x[n-2], "; ", zeta, "=", v), list(v=lags[3]))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.45,0.5), breaks=c(seq(-0.4,-0.1,0.1), 0, seq(0.1,0.5,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot8_3


#### zeta = 1

### 1) results for xt
res <- reshape(df3_6, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16,17:19,20:22,23:25,26:28,29:31))
res$up <- res$bias1 + 2*res$biassd1 
res$lo <- res$bias1 - 2*res$biassd1 

res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot8_4 <- ggplot(res, aes(phi1, bias1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias in % of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste("Bias in ", x[n], "; ", zeta, "=", v), list(v=lags[4]))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.45,0.5), breaks=c(seq(-0.4,-0.1,0.1), 0, seq(0.1,0.5,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot8_4


### 2) results for xt-1
res <- reshape(df3_6, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16,17:19,20:22,23:25,26:28,29:31))
res$up <- res$bias1l1 + 2*res$biassd1l1 
res$lo <- res$bias1l1 - 2*res$biassd1l1 

res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot8_5 <- ggplot(res, aes(phi1, bias1l1, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias in % of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste("Bias in ", x[n-1], "; ", zeta, "=", v), list(v=lags[4]))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.45,0.5), breaks=c(seq(-0.4,-0.1,0.1), 0, seq(0.1,0.5,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot8_5


### 3) results for xt-2
res <- reshape(df3_6, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16,17:19,20:22,23:25,26:28,29:31))
res$up <- res$bias1l2 + 2*res$biassd1l2 
res$lo <- res$bias1l2 - 2*res$biassd1l2 

res$time <- factor(res$time, levels=c(1,2,3), labels=c("RE model","FE model","FEIS model"))
#View(res)

simplot8_6 <- ggplot(res, aes(phi1, bias1l2, group=time, color=time, shape=time)) +
  geom_line() +
  geom_pointrange(aes(ymin = lo, ymax = up, color=time)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #  geom_vline(xintercept = 0, color="red") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab(expression(paste("Mean Bias in % of ",beta,"  +/- 2 SD"))) +
  ggtitle(substitute(paste("Bias in ", x[n-2], "; ", zeta, "=", v), list(v=lags[4]))) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  scale_y_continuous(limits=c(-0.45,0.5), breaks=c(seq(-0.4,-0.1,0.1), 0, seq(0.1,0.5,0.1))) +
  theme_bw() +
  theme(legend.key=element_blank(),
        legend.position=c(0.1,0.9), legend.justification=c(0.1,0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot8_6



### Combine Figures
library(ggpubr)
simplot8 <- ggarrange(simplot8_1,simplot8_2,simplot8_3,simplot8_4,simplot8_5,simplot8_6,
                      ncol=3, nrow=2, common.legend = F, 
                      labels=c("A) ", "B) ", "C) ", "D) ", "E) ", "F) "), font.label = list(size = 20)) 

cairo_pdf(file=paste("../03_Output/Output_", "simplot7_lag2.pdf", sep=""), width=17.54, height=10, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
simplot8
dev.off()





#### Rejection rate

### 1) results for zeta = 0.1
res <- reshape(df0_6, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16,17:19,20:22,23:25,26:28,29:31))

res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot9_1 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste(zeta, "=", v), list(v=lags[1]))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot9_1

### 2) results for zeta = 0.3
res <- reshape(df1_6, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16,17:19,20:22,23:25,26:28,29:31))

res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot9_2 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste(zeta, "=", v), list(v=lags[2]))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot9_2


### 3) results for zeta = 0.5
res <- reshape(df2_6, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16,17:19,20:22,23:25,26:28,29:31))

res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot9_3 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste(zeta, "=", v), list(v=lags[3]))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot9_3


### 4) results for zeta = 1
res <- reshape(df3_6, direction="long", varying=list(2:4,5:7,8:10,11:13,14:16,17:19,20:22,23:25,26:28,29:31))

res$time <- factor(res$time, levels=c(1,2,3), labels=c("FEIS vs. FE","FE vs RE","FEIS vs. RE"))

simplot9_4 <- ggplot(res, aes(phi1,rej1, group=time, color=time, shape=time)) +
  geom_point(size=3) +
  geom_line() +
  geom_hline(yintercept = 0.05,color="red", linetype=2) +
  geom_hline(yintercept = 0.95,color="red", linetype=2) +
  geom_vline(xintercept = 0,color="black") +
  xlab(expression(paste("Cov(",delta," , ",alpha[2],")","     (",phi,")"))) +
  ylab("Proportion rejected") +
  ggtitle(substitute(paste(zeta, "=", v), list(v=lags[4]))) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(-0.85,0.85), 
                     breaks = c(seq(-0.8,0.8,0.2))) +
  theme_bw()  +
  theme(legend.key=element_blank(),
        legend.position=c(0.95,0.4), legend.justification=c(0.95,0.4),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, margin = margin(r=10)),
        plot.title = element_text(size=20),
        strip.background =element_blank(),
        strip.text = element_text(size=20, colour="black"),
        panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())

simplot9_4



### Combine Figures
library(ggpubr)
simplot9 <- ggarrange(simplot9_1,simplot9_2,simplot9_3,simplot9_4,
                      ncol=2, nrow=2, common.legend = F, 
                      labels=c("A) ", "B) ", "C) ", "D) "), font.label = list(size = 20)) 

cairo_pdf(file=paste("../03_Output/Output_", "simplot7_lag2_rejection.pdf", sep=""), width=17.54, height=10, 
          bg = "white", family="Times New Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
simplot9
dev.off()
