#### ---------------------------------------------------------------####
# Example 1 : Male Marital wage premium 
# R Script for article 
# Rüttenauer T. and Ludwig V. 2020. Fixed Effects Individual Slopes: 
# Accounting and Testing for Heterogeneous Effects in Panel Data 
# or other Multilevel Models. Sociological Methods and Research
#### ---------------------------------------------------------------####



####################
#### Example 1 #####
####################


rm(list=ls())


#### load package

library(feisr)
library(plm)
library(lme4)
library(haven)
library(texreg)

#### Set working directory
setwd("")



#### Load data
# This dataset is comes for Ludwig / Brüderl replication files
mwp.df <- data.frame(haven::read_dta("../01_Stata/mwp_US_analysis.dta"))



######################
#### Prepare data ####
######################

### Recode child and year 
mwp.df$dchild <- mwp.df$nchild
mwp.df$dchild[mwp.df$dchild == 0] <- 0
mwp.df$dchild[mwp.df$dchild == 1] <- 1
mwp.df$dchild[mwp.df$dchild == 2] <- 2
mwp.df$dchild[mwp.df$dchild >= 3] <- 3

mwp.df$yeargr <- mwp.df$year
mwp.df$yeargr[mwp.df$yeargr %in% c(1979:1980)] <- 1
mwp.df$yeargr[mwp.df$yeargr %in% c(1981:1985)] <- 2
mwp.df$yeargr[mwp.df$yeargr %in% c(1986:1990)] <- 3
mwp.df$yeargr[mwp.df$yeargr %in% c(1991:1995)] <- 4
mwp.df$yeargr[mwp.df$yeargr %in% c(1996:2000)] <- 5
mwp.df$yeargr[mwp.df$yeargr %in% c(2001:2005)] <- 6
mwp.df$yeargr[mwp.df$yeargr %in% c(2006:2012)] <- 7




###########################
#### Regression models ####
###########################


### Select sample
mwp.df <- mwp.df[mwp.df$sample1==1, ]


### Estimate FEIS non-robust SEs
wages.feis <- feis(lnw ~ marry + enrol + yeduc + as.factor(dchild) + tenure + as.factor(yeargr)
                       | exp + I(exp^2), data = mwp.df, id = "id",
                       robust = FALSE)
summary(wages.feis)

### Estimate FEIS robust
wages.feis.rob <- feis(lnw ~ marry + enrol + yeduc + as.factor(dchild) + tenure + as.factor(yeargr)
                       | exp + I(exp^2), data = mwp.df, id = "id",
                       robust = TRUE)
summary(wages.feis.rob)



### Estimate FE and RE
wages.fe <- plm(lnw ~ marry + enrol + yeduc + as.factor(dchild) + tenure + as.factor(yeargr)
                + exp + I(exp^2), data = mwp.df, index = c("id", "year"),
                model = "within", effect = "individual")
wages.re <- plm(lnw ~ marry + enrol + yeduc + as.factor(dchild) + tenure + as.factor(yeargr)
                + exp + I(exp^2), data = mwp.df, index = c("id", "year"),
                model = "random", effect = "individual")

# Compute robust SEs
wages.fe.rob <- wages.fe
wages.fe.rob$vcov <- vcovHC(wages.fe.rob, method = "arellano", cluster="group") 
summary(wages.fe.rob)

wages.re.rob <- wages.re
wages.re.rob$vcov <- vcovHC(wages.re.rob, method = "arellano", cluster="group")
summary(wages.re.rob)



### Estimate RE RS model
wages.rs <- lmer(lnw ~ marry + enrol + yeduc + as.factor(dchild) + tenure + as.factor(yeargr) + exp + I(exp^2) +
             (1 + exp + I(exp^2) | id), data = mwp.df)
summary(wages.rs)




########################
#### Output Table 1 ####
########################


# Set method for texreg
setMethod("extract", signature = className("feis", "feisr"),
          definition = feisr:::extract.feis)

screenreg(list(wages.re.rob, wages.fe.rob, wages.rs, wages.feis.rob), digits = 3,
          custom.model.names = c("RE", "FE", "RERS", "FEIS"))


### Specify file name
file <- "Table_Example1.tex"

# Save regression table
tab <- texreg(l = list(wages.re.rob, wages.fe.rob, wages.rs, wages.feis.rob), 
              #file = file,
              digits = 3, leading.zero = TRUE,
              stars = c(0.001, 0.01, 0.05),
              symbol = "\\cdot",
              caption = "Example 1. Dep var: ln hourly wage rate",
              custom.coef.names = c("Married", "Currently enrolled", "Years education", 
                                    "One child", "Two children", "Three or more children", "Tenure (years)",
                                     # "Age group 2", "Age group 3", "Age group 4", "Age group 5", "Age group 6", "Age group 7",
                                    "Work experience", "Work experience$^2$"),
              custom.model.names = c("RE", "FE", "RIRS", "FEIS"),
              omit.coef = "(Intercept)|(yeargr)", label = "table:example1",
              custom.note = "%stars. Robust standard errors in parentheses. Age group omitted.",
              dcolumn = TRUE, caption.above = TRUE, use.packages = FALSE, float.pos = "t",
              include.aic = FALSE, include.bic = FALSE,
              include.loglik = FALSE, include.variance = FALSE, include.rmse = FALSE)

# Customize column space
tab <- gsub("D[{].[}][{].[}][{][[:digit:]]+\\.*[[:digit:]]*[}]", "D{.}{.}{4.6}", tab)

# Center N and god
l1 <- gregexpr("R$^2$", tab, fixed = TRUE) #start gofs
l1 <- l1[[1]][1]

l2 <- gregexpr("\\hline", tab, fixed = TRUE) #end gofs
l2 <- l2[[1]][length(l2[[1]])]

tmp <- substr(tab, l1, (l2-1))
tmp2 <- gsub("& ([[:digit:]]+|[[:digit:]]+\\.[[:digit:]]+)", "& \\\\multicolumn{1}{c}{\\1}", tmp)
tab <- gsub(tmp, tmp2, tab, fixed = TRUE)

# Export
write.table(tab, file = file,            
            col.names = FALSE, row.names = FALSE, quote = FALSE)




#####################################
#### Artificial regression tests ####
#####################################


### Artificial regression test
ht <- feistest(wages.feis.rob, robust = TRUE, type = "all")
summary(ht)


### Bootstrapped Hausman test
bsht <- bsfeistest(wages.feis.rob, type = "all", rep = 100, seed = 123)
summary(bsht)




########################
#### Output Table 2 ####
########################

### Specify file name
file <- "Spec_Example1.tex"


### Write latex table
cat("\\begin{table}\n", file = file)
cat("\\caption{Example 1. Specification tests}\n", file = file, append = TRUE)
cat("\\label{table:example1_spec}\n", file = file, append = TRUE)
cat("\\centering\n", file = file, append = TRUE)
cat("\\begin{tabular}{l D{.}{.}{2.3} D{.}{.}{2.0} D{.}{.}{2.6} }", file = file, append = TRUE)
cat("\\hline \n", file = file, append = TRUE)
cat(" & \\multicolumn{1}{c}{$\\chi^2$} & \\multicolumn{1}{c}{df} & \\multicolumn{1}{c}{P(> $\\chi^2$)}  \\\\ \n", file = file, append = TRUE)
cat("\\hline \n", file = file, append = TRUE)

# ART results
cat("Artificial regression test \\\\ \n", file = file, append = TRUE)

cat("FEIS vs. FE:", 
    format(unname(ht$wald_feis$result$chi2[1]), digits = 3, nsmall = 3), 
    format(unname(ht$wald_feis$result$chi2[2]), digits = 3),
    format(unname(ht$wald_feis$result$chi2[3]), digits = 3, nsmall = 3),
    sep = " & "
    , file = file, append = TRUE) 
cat(" \\\\ \n", file = file, append = TRUE) 

cat("FE vs. RE:", 
    format(unname(ht$wald_fe$result$chi2[1]), digits = 3, nsmall = 3),  
    format(unname(ht$wald_fe$result$chi2[2]), digits = 3), 
    format(unname(ht$wald_fe$result$chi2[3]), digits = 3, nsmall = 3),  
    sep = " & "
    , file = file, append = TRUE) 
cat(" \\\\ \n", file = file, append = TRUE) 

cat("FEIS vs. RE:", 
    format(unname(ht$wald_re$result$chi2[1]), digits = 3, nsmall = 3),  
    format(unname(ht$wald_re$result$chi2[2]), digits = 3), 
    format(unname(ht$wald_re$result$chi2[3]), digits = 3, nsmall = 3),  
    sep = " & "
    , file = file, append = TRUE) 
cat(" \\\\ \n", file = file, append = TRUE) 


# BSH results
cat(" \\\\ \n", file = file, append = TRUE)
cat("Bootstrapped Hausman test \\\\ \n", file = file, append = TRUE)

cat("FEIS vs. FE:", 
    format(unname(bsht$wald_feis$result$chi2[1]), digits = 3, nsmall = 3), 
    format(unname(bsht$wald_feis$result$chi2[2]), digits = 3),
    format(unname(bsht$wald_feis$result$chi2[3]), digits = 3, nsmall = 3),
    sep = " & "
    , file = file, append = TRUE) 
cat(" \\\\ \n", file = file, append = TRUE) 

cat("FE vs. RE:", 
    format(unname(bsht$wald_fe$result$chi2[1]), digits = 3, nsmall = 3),  
    format(unname(bsht$wald_fe$result$chi2[2]), digits = 3), 
    format(unname(bsht$wald_fe$result$chi2[3]), digits = 3, nsmall = 3),  
    sep = " & "
    , file = file, append = TRUE) 
cat(" \\\\ \n", file = file, append = TRUE) 

cat("FEIS vs. RE:", 
    format(unname(bsht$wald_re$result$chi2[1]), digits = 3, nsmall = 3),  
    format(unname(bsht$wald_re$result$chi2[2]), digits = 3), 
    format(unname(bsht$wald_re$result$chi2[3]), digits = 3, nsmall = 3),  
    sep = " & "
    , file = file, append = TRUE) 
cat(" \\\\ \n", file = file, append = TRUE) 

# End
cat("\\hline\n", file = file, append = TRUE)
cat("\\end{tabular}\n", file = file, append = TRUE)
cat("\\end{table}\n", file = file, append = TRUE)



