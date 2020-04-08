
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduction-Material-Fixed-Effects-Individual-Slopes

This repository provides all materials for reproduction of results in

> Rüttenauer T. and Ludwig V. (2020). “Fixed Effects Individual Slopes:
> Accounting and Testing for Heterogeneous Effects in Panel Data or
> other Multilevel Models”, *Sociological Methods and Research*,
> Forthcoming.

## Requirements

The codes needs the following folders: “01\_Script”, “02\_Data”,
“03\_Output”. All R Scripts are required in folder “01\_Script”.

The following packages are necessary for reproduction of main results:

``` r
install.packages("Formula")
install.packages("Matrix")
install.packages("MASS")
install.packages("matrixcalc")
install.packages("plm")
install.packages("aod")
install.packages("feisr") 
install.packages("doParallel")
install.packages("ggplot2")
install.packages("lattice")
install.packages("gridExtra")
install.packages("ggpubr")
install.packages("extrafont")
```

## Replication: Monte Carlo simulations

The R script named “00\_FEIS\_Simulation\_Program.R” contains the
simulation program (function *fesim*). The R script
“00\_FEIS\_Simulation\_Program\_hpc.R” contains the equivalent
function for parallelization on a high performance cluster. The R script
“01\_FEIS\_Simulation\_Main-Analysis.R” replicates the results of the
main paper. To run the replication files on a HPC, comment line 85
(loading local version of the program) and uncomment line 88 (loading
the HPC version of the program). Furthermore, comment the working (line
57) directory, if specified in an external batch file. The script
“02\_FEIS\_Simulation\_Supplementary-Analyses.R” replicates the
supplementary results. The batch files “01\_runR.sl” and “02\_runR.sl”
can be used to run the analyses on a Slurm HPC.

The folder “02\_Data” contains the results of individual simulation runs
as performed and saved in the replication files.

Parameters (number of observations, periods, replication and bootstrap
replications) can be set at the beginning of the replication files. Note
that the supplementary analyses define to different sets of parameters,
once around line 85 for small N simulations, once around line 1287 for
lagged treatment effects and panel selection.

## Replication: example 1

tbd

Replication material of Ludwig, V., & Brüderl, J. (2018) "Is There a Male Marital Wage Premium? New Evidence from the United States" cam be obtained from:
<https://github.com/volkerludwig/Replication-material-for-Ludwig-Bruederl-ASR2018>

## Replication: example 2

tbd

Replication material of Deming, D. (2009) “Early Childhood Intervention
and Life-Cycle Skill Development: Evidence from Head Start” can be
obtained from:
<https://www.aeaweb.org/maintenance.php?id=10.1257/app.1.3.111>

## Run time of simulations

Note that using 1000 replications (reps = 1000) and 100 bootstrap
replications (bsreps = 100) takes over 40 hours on 96 parallel cores.
For replication purposes, we recommend to reduce the number of
simulations trials or bootstrap runs. Especially increasing the number
of bootstrap runs in the BHT tests is computationally intense.

## System and version information

Platform: Windows 10 (x86\_64) | GNU/Linux SMP (x86\_64)

Version: R version 3.6.1
