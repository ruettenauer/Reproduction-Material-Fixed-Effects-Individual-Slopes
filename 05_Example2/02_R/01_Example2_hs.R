#### ---------------------------------------------------------------####
# Example 2 : Effectiveness of the Head Start preschool program
# R Script for article 
# Rüttenauer T. and Ludwig V. 2020. Fixed Effects Individual Slopes: 
# Accounting and Testing for Heterogeneous Effects in Panel Data 
# or other Multilevel Models. Sociological Methods and Research
#### ---------------------------------------------------------------####



####################
#### Example 2 #####
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
# This dataset is created by the Stata dofile "01_Example2_hs.do" based on Demings replication files
nlsy.df <- data.frame(haven::read_dta("../01_Stata/hs.dta"))



###########################
#### Regression models ####
###########################

# Gen fake time
nlsy.df$ft <- ave(1:nrow(nlsy.df), by = nlsy.df$MotherID, FUN = function(x) 1:length(x))



# Replicate M5 Table 3 from Deming
rep.fe <- plm(Test_std ~ HS_5to6 + HS_7to10 + HS_11to14 + Pre_5to6 + Pre_7to10 + Pre_11to14 +
              Male + as.factor(year) + Group_5to6 + Group_7to10 + Group_11to14 + as.factor(AgeTest_Yr) +
              Res_0to3_imp + VLow_BW_imp + logBW_imp + LogInc_0to3_imp + LogIncAt3_imp + FirstBorn_imp + 
              PPVTat3_imp + HOME_Pct_0to3_imp + Moth_HrsWorked_BefBirth_imp + Moth_HrsWorked_Avg_0to3_imp + 
              Moth_HrsWorked_0to1_imp + Father_HH_0to3_imp + MomCare_imp + NonRelCare_imp + Breastfed_imp + 
              Doctor_0to3_imp + Dentist_0to3_imp + Moth_WeightChange_imp + Insurance_0to3_imp +
              Res_0to3_miss + VLow_BW_miss + logBW_miss + LogInc_0to3_miss + LogIncAt3_miss + FirstBorn_miss + 
              PPVTat3_miss + HOME_Pct_0to3_miss + Moth_HrsWorked_BefBirth_miss + Moth_HrsWorked_Avg_0to3_miss + 
              Moth_HrsWorked_0to1_miss + Father_HH_0to3_miss + MomCare_miss + NonRelCare_miss + Breastfed_miss + 
              Doctor_0to3_miss + Dentist_0to3_miss + Moth_WeightChange_miss + Insurance_0to3_miss +
              HealthCond_before_imp + GMom_0to3_imp + MomCare_imp + RelCare_imp + Moth_Smoke_BefBirth_imp + 
              Alc_BefBirth_imp + Illness_1stYr_imp + Premature_imp + Medicaid_0to3_imp +
              HealthCond_before_miss + GMom_0to3_miss + MomCare_miss + RelCare_miss + Moth_Smoke_BefBirth_miss + 
              Alc_BefBirth_miss + Illness_1stYr_miss + Premature_miss + Medicaid_0to3_miss
                , 
                data = nlsy.df, model = "within", effect = "individual",
                index = c("MotherID", "ft"))
rep.fe$vcov <- vcovHC(rep.fe, method = "arellano", cluster = "group")
summary(rep.fe)


# Get FEIS sample 
sample <- feis(Test_std ~ HS_5to6 + HS_7to10 + HS_11to14 + Pre_5to6 + Pre_7to10 + Pre_11to14 +
                 as.factor(year) + Group_5to6 + Group_7to10 + Group_11to14 + as.factor(AgeTest_Yr)
               | PreTreatIndex, 
               data = nlsy.df, id = "MotherID", robust = T)
nlsy.df$sample <- 1
nlsy.df$sample[sample$na.action] <- 0


# Replicate M5 with Pretreatment Index
head1.fe <- plm(Test_std ~ HS_5to6 + HS_7to10 + HS_11to14 + Pre_5to6 + Pre_7to10 + Pre_11to14 +
                as.factor(year) + Group_5to6 + Group_7to10 + Group_11to14 + as.factor(AgeTest_Yr) +
                PreTreatIndex
               , 
               data = nlsy.df[nlsy.df$sample == 1, ], 
               model = "within", effect = "individual",
               index = c("MotherID", "ft"))
head1.fe$vcov <- vcovHC(head1.fe, method = "arellano", cluster = "group")
summary(head1.fe)


# FEIS model of head1
head1.feis <- feis(Test_std ~ HS_5to6 + HS_7to10 + HS_11to14 + Pre_5to6 + Pre_7to10 + Pre_11to14 +
                    as.factor(year) + Group_5to6 + Group_7to10 + Group_11to14 + as.factor(AgeTest_Yr)
                  | PreTreatIndex, 
                  data = nlsy.df, id = "MotherID", robust = T)
summary(head1.feis)


#### Reduce dataset to families with Head Start Variance

nlsy2.df <- nlsy.df[nlsy.df$sample==1, ]

nlsy2.df$nchild <- ave(nlsy2.df$ChildID, nlsy2.df$MotherID, FUN = function(x) length(unique(x)))

nlsy2.df$sd <- ave(nlsy2.df$HS2_FE90, nlsy2.df$MotherID, FUN = function(x) sd(x))

nlsy2.df$nobs <- ave(nlsy2.df$HS2_FE90, nlsy2.df$MotherID, FUN = function(x) length(x))

nlsy2.df <- nlsy2.df[which(nlsy2.df$nchild >= 2 & nlsy2.df$sd != 0), ] 

# Fake year
nlsy2.df$ft <- ave(1:nrow(nlsy2.df), by = nlsy2.df$MotherID, FUN = function(x) 1:length(x))


# Get FEIS sample 
sample2 <- feis(Test_std ~ HS_5to6 + HS_7to10 + HS_11to14 + Pre_5to6 + Pre_7to10 + Pre_11to14 +
                 as.factor(year) + Group_5to6 + Group_7to10 + Group_11to14 + as.factor(AgeTest_Yr)
               | PreTreatIndex, 
               data = nlsy2.df, id = "MotherID", robust = T)
nlsy2.df$sample2 <- 1
nlsy2.df$sample2[sample2$na.action] <- 0

# Replicate M5 with Pretreatment Index
head2.fe <- plm(Test_std ~ HS_5to6 + HS_7to10 + HS_11to14 + Pre_5to6 + Pre_7to10 + Pre_11to14 +
                  as.factor(year) + Group_5to6 + Group_7to10 + Group_11to14 + as.factor(AgeTest_Yr) 
                + PreTreatIndex
                , 
                data=nlsy2.df[nlsy2.df$sample2 == 1, ], 
                model = "within", effect = "individual",
                index = c("MotherID", "ft"))
head2.fe$vcov <- vcovHC(head2.fe, method = "arellano", cluster = "group")
summary(head2.fe)


# FEIS model of head2
head2.feis <- feis(Test_std ~ HS_5to6 + HS_7to10 + HS_11to14 + Pre_5to6 + Pre_7to10 + Pre_11to14 +
                     as.factor(year) + Group_5to6 + Group_7to10 + Group_11to14 + as.factor(AgeTest_Yr)
                   | PreTreatIndex, 
                   data = nlsy2.df, id = "MotherID", robust = T)
summary(head2.feis)

# Random Slopes model of head2
head2.rs <- lmer(Test_std ~ HS_5to6 + HS_7to10 + HS_11to14 + Pre_5to6 + Pre_7to10 + Pre_11to14 +
                   as.factor(year) + Group_5to6 + Group_7to10 + Group_11to14 + as.factor(AgeTest_Yr) +
                   PreTreatIndex +
                   (1 + PreTreatIndex | MotherID), data = nlsy2.df)
summary(head2.rs)




########################
#### Output Table 3 ####
########################

# Set method for texreg
setMethod("extract", signature = className("feis", "feisr"),
          definition = feisr:::extract.feis)

screenreg(list(rep.fe, head1.fe, head1.feis, head2.fe, head2.rs, head2.feis), digits = 3,
          custom.model.names = c("Rep", "FE", "FEIS", "FE 2", "RERS 2", "FEIS 2"))

### Specify file name
file <- "Table_Example2.tex"

### Save latex table
tab <- texreg(l = list(rep.fe, head1.fe, head1.feis, head2.fe, head2.rs, head2.feis), 
              # file = file,
              digits = 3, leading.zero = TRUE,
              stars = c(0.001, 0.01, 0.05),
              symbol = "\\cdot",
              caption = "Example 2. Dep var: cognitive test scores",
              custom.coef.names = c("Ages 5-6", "Ages 7-10", "Ages 11-14",  
                                    "Ages 5-6 ", "Ages 7-10 ", "Ages 11-14 ", 
                                    "Pretreatment index"),
              groups = list("Head Start" = 1:3, "Other Preschool" = 4:6),
              custom.model.names = c("Rep", "FE", "FEIS", "FE 2", "RIRS 2", "FEIS 2"),
              omit.coef = "(year)|(AgeTest_Yr)|(Group_)|(_imp)|(_miss)|(Intercept)|(Male)",
              label = "table:example2",
              custom.note = "%stars. Robust standard errors in parentheses.",
              dcolumn = TRUE, caption.above = TRUE, use.packages = FALSE, float.pos = "t",
              include.aic = FALSE, include.bic = FALSE,
              include.loglik = FALSE, include.variance = FALSE, include.rmse = FALSE)

# Customize column space
tab <- gsub("D[{].[}][{].[}][{][[:digit:]]+\\.*[[:digit:]]*[}]", "D{.}{.}{2.4}", tab)

# Customize header
oh <- "& \\multicolumn{1}{c}{Rep} & \\multicolumn{1}{c}{FE} & \\multicolumn{1}{c}{FEIS} & \\multicolumn{1}{c}{FE 2} & \\multicolumn{1}{c}{RIRS 2} & \\multicolumn{1}{c}{FEIS 2} \\\\"
nh <- paste0("& \\multicolumn{3}{c}{Orig. Sample} & \\multicolumn{3}{c}{Subsample} \\\\ \n",
             "\\cmidrule(r){2-4} \\cmidrule(l){5-7} \n",
             oh, "\n",
             "& \\multicolumn{1}{c}{(1)} & \\multicolumn{1}{c}{(2)} & \\multicolumn{1}{c}{(3)} & \\multicolumn{1}{c}{(4)} & \\multicolumn{1}{c}{(5)} & \\multicolumn{1}{c}{(6)} \\\\")
tab <- gsub(oh, nh, tab, fixed = TRUE)

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


### Artificial regression test original sample
ht1 <- feistest(head1.feis, robust = TRUE, type = "all")
summary(ht1)


### Bootstrapped Hausman test original sample
bsht1 <- bsfeistest(head1.feis, type = "all", rep = 100, seed = 123)
summary(bsht1)


### Artificial regression test subsample
ht <- feistest(head2.feis, robust = TRUE, type = "all")
summary(ht)


### Bootstrapped Hausman test subsample
bsht <- bsfeistest(head2.feis, type = "all", rep = 100, seed = 123)
summary(bsht)



########################
#### Output table 4 ####
########################

### Specify file name
file <- "Spec_Example2.tex"


### Write latex table
cat("\\begin{table}\n", file = file)
cat("\\caption{Example 2. Specification tests}\n", file = file, append = TRUE)
cat("\\label{table:example2_spec}\n", file = file, append = TRUE)
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
    format(unname(ht$wald_feis$result$chi2[3]), digits = 2, nsmall = 3),
    sep = " & "
    , file = file, append = TRUE) 
cat(" \\\\ \n", file = file, append = TRUE) 

cat("FE vs. RE:", 
    format(unname(ht$wald_fe$result$chi2[1]), digits = 3, nsmall = 3),  
    format(unname(ht$wald_fe$result$chi2[2]), digits = 3), 
    format(unname(ht$wald_fe$result$chi2[3]), digits = 2, nsmall = 3),  
    sep = " & "
    , file = file, append = TRUE) 
cat(" \\\\ \n", file = file, append = TRUE) 

cat("FEIS vs. RE:", 
    format(unname(ht$wald_re$result$chi2[1]), digits = 3, nsmall = 3),  
    format(unname(ht$wald_re$result$chi2[2]), digits = 3), 
    format(unname(ht$wald_re$result$chi2[3]), digits = 2, nsmall = 3),  
    sep = " & "
    , file = file, append = TRUE) 
cat(" \\\\ \n", file = file, append = TRUE) 


# BSH results
cat(" \\\\ \n", file = file, append = TRUE)
cat("Bootstrapped Hausman test \\\\ \n", file = file, append = TRUE)

cat("FEIS vs. FE:", 
    format(unname(bsht$wald_feis$result$chi2[1]), digits = 3, nsmall = 3), 
    format(unname(bsht$wald_feis$result$chi2[2]), digits = 3),
    format(unname(bsht$wald_feis$result$chi2[3]), digits = 1, nsmall = 1),
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
    format(unname(bsht$wald_re$result$chi2[3]), digits = 2, nsmall = 2),  
    sep = " & "
    , file = file, append = TRUE) 
cat(" \\\\ \n", file = file, append = TRUE) 

# End
cat("\\hline\n", file = file, append = TRUE)
cat("\\end{tabular}\n", file = file, append = TRUE)
cat("\\end{table}\n", file = file, append = TRUE)














