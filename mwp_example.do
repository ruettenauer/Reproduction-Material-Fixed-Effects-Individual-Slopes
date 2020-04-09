*=======================================================================
* Example 1 : Male Marital wage premium 
* Stata Do-file for article 
* Rüttenauer T. and Ludwig V. 2020. Fixed Effects Individual Slopes: 
* Accounting and Testing for Heterogeneous Effects in Panel Data 
* or other Multilevel Models. Sociological Methods and Research
*=======================================================================

* This Stata 15 Do-file replicates the analysis by Ludwig and Brüderl (2018)

//set working directory
global xmpldir ""
cd $xmpldir


*=======================================================================
* Example
*=======================================================================

use mwp_US_analysis, clear
xtset id year
xtsum lnw exp marry evermarr nchild enrol yeduc year

recode nchild 3/max=3, gen(dchild)
tab dchild, gen(dchild)
recode year 1979/1980=1 1981/1985=2 1986/1990=3 1991/1995=4 1996/2000=5 2001/2005=6 2006/2012=7, gen(yeargr)
tab yeargr, gen(yeargr)
	
keep if sample1==1

//FEIS estimation and comparison across models

xtfeis lnw marry dchild2-dchild4 enrol yeduc tenure ib1.yeargr, slope(c.exp c.exp#c.exp) cluster(id)
est store FEIS

//standard FE model
qui xtreg lnw marry dchild2-dchild4  enrol yeduc tenure c.exp##c.exp ib1.yeargr, fe vce(cluster id) 
qui est store FE

//reproduce FE model using xtfeis command
qui xtfeis lnw marry dchild2-dchild4  enrol yeduc tenure c.exp##c.exp ib1.yeargr, cluster(id) 
qui est store FE_xtfeis

//FEGS model (FE allowing for Group-specific Slopes)
qui xtreg lnw marry dchild2-dchild4  enrol yeduc tenure c.exp c.exp#c.exp ib1.yeargr ///
	1.evermarr#c.exp 1.evermarr#c.exp#c.exp , fe vce(cluster id)
qui est store FEGS

//Standard RE model
qui xtreg lnw marry dchild2-dchild4  enrol yeduc tenure c.exp##c.exp ib1.yeargr , re vce(cluster id) 
qui est store RE

//RIRS model
qui mixed lnw marry dchild2-dchild4  enrol yeduc tenure ib1.yeargr exp expq || id: exp expq , vce(cluster id) mle
qui est store RS

*qui mixed lnw marry dchild2-dchild4  enrol yeduc tenure ib1.yeargr exp expq || id: exp , vce(cluster id) mle
*qui est store RS

//RSGS model
qui mixed lnw marry dchild2-dchild4  enrol yeduc tenure ib1.yeargr exp expq evermarr 1.evermarr#c.exp 1.evermarr#c.expq || id: exp expq , vce(cluster id) 
qui est store RSGS
test evermarr
test 1.evermarr#c.exp 1.evermarr#c.expq 


estout RE FE FEIS RS FEGS RSGS, cells(b(fmt(3)) se(par fmt(3))) ///
	numbers varwidth(22) stats(N N_g r2_w, fmt(0 0 3) ///
	labels("Number of person-years" "Number of persons" "R^2 within")) ///
	keep(marry dchild2 dchild3 dchild4 enrol yeduc tenure exp c.exp#c.exp 1.evermarr#c.exp ///
	1.evermarr#c.exp#c.exp expq 1.evermarr#c.expq ) ///
	note("Note: Coefficients and clustered standard errors in parentheses.") 

file open Spec using Table_Example1.tex, replace write

file write Spec "\begin{table}[t]" _n
file write Spec "\caption{Example 1: The marital wage premium for men}" _n
file write Spec "\begin{center}" _n
file write Spec "\begin{tabular}{l D{.}{.}{4.6} D{.}{.}{4.6} D{.}{.}{4.6} D{.}{.}{4.6} }" _n
file write Spec "\hline & \multicolumn{1}{c}{RE} & \multicolumn{1}{c}{FE} & \multicolumn{1}{c}{FEIS} & \multicolumn{1}{c}{RIRS} \\" _n
file write Spec "\hline" _n

file close Spec

	
estout RE FE FEIS RS using Table_Example1.tex, cells(b(fmt(3)) se(par fmt(3))) ///
	numbers varwidth(22) stats(N N_g r2_w, fmt(0 0 3) ///
	labels("Number of person-years" "Number of persons" "R$^2$ within")) ///
	keep(marry dchild2 dchild3 dchild4 enrol yeduc tenure exp c.exp#c.exp) ///
	varlabels(_cons \_cons)	///
	note("Note: Coefficients and cluster-robust standard errors in parentheses.") ///
	style(tex) append


	
// Artificial Regression Test (ART) and Bootstrapped Hausman Test (BSHT) 

//(xtart requires FEIS model specified w/o using factor notation)
qui xtfeis lnw marry dchild2-dchild4 enrol yeduc tenure yeargr2-yeargr7, slope(exp expq) cluster(id)
est store FEIS

//standard FE model
qui xtreg lnw marry dchild2-dchild4  enrol yeduc tenure exp expq yeargr2-yeargr7, fe vce(cluster id) 
qui est store FE

//Standard RE model
qui xtreg lnw marry dchild2-dchild4  enrol yeduc tenure exp expq yeargr2-yeargr7, re vce(cluster id) 
qui est store RE


xtart FEIS
di "chi2(`r(df)') = `r(chi2)'"
di "p-value = `r(p)'"


xtart FEIS
local ARTchi2=round(`r(chi2)',0.001)
local ARTdf=`r(df)'
local ARTp=`r(p)' 
xtart FEIS, fe
local ARTchi2_1=round(`r(chi2)',0.001) 
local ARTdf_1=`r(df)' 
local ARTp_1=`r(p)' 
xtart FEIS, re
local ARTchi2_2=round(`r(chi2)',0.001) 
local ARTdf_2=`r(df)' 
local ARTp_2=`r(p)' 

xtbsht FEIS FE, seed(123) reps(100)
local BSHTchi2=round(`r(chi2)',0.001)
local BSHTdf=`r(df)'
local BSHTp=`r(p)' 
xtbsht FE RE, seed(123) reps(100)
local BSHTchi2_1=round(`r(chi2)',0.001) 
local BSHTdf_1=`r(df)' 
local BSHTp_1=`r(p)' 
xtbsht FEIS RE, seed(123) reps(100)
local BSHTchi2_2=round(`r(chi2)',0.001) 
local BSHTdf_2=`r(df)' 
local BSHTp_2=`r(p)' 

file open Spec using Spec_Example1.tex, replace write

file write Spec "\begin{table}" _n
file write Spec "\caption{Example 1: Specification tests}" _n
file write Spec "\label{table:example1_spec}" _n
file write Spec "\centering" _n
file write Spec "\begin{tabular}{l D{.}{.}{2.3} D{.}{.}{2.0} D{.}{.}{2.6} }\hline  & \multicolumn{1}{c}{$\chi^2$} & \multicolumn{1}{c}{df} & \multicolumn{1}{c}{P(> $\chi^2$)}  \\" _n 
file write Spec "\hline" _n 
file write Spec "Artificial regression test \\" _n 
file write Spec "FEIS vs. FE: & `ARTchi2' & `ARTdf' & `ARTp' \\" _n 
file write Spec "FE vs. RE: & `ARTchi2_1' & `ARTdf_1' & `ARTp_1' \\" _n 
file write Spec "FEIS vs. RE: & `ARTchi2_2' & `ARTdf_2' & `ARTp_2' \\" _n 
file write Spec " \\" _n 
file write Spec "Bootstrapped Hausman test \\" _n 
file write Spec "FEIS vs. FE: & `BSHTchi2' & `BSHTdf' & `BSHTp' \\" _n 
file write Spec "FE vs. RE: & `BSHTchi2_1' & `BSHTdf_1' & `BSHTp_1' \\" _n 
file write Spec "FEIS vs. RE: & `BSHTchi2_2' & `BSHTdf_2' & `BSHTp_2' \\" _n 
file write Spec "\hline" _n
file write Spec "\end{tabular}" _n
file write Spec "\end{table}" _n

file close Spec



