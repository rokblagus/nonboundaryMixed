This is all the source code needed to reproduce the results repoted in the paper "Non-boundary covariance matrix estimation in generalized linear mixed effects models using data augmentation priors", by Tina Košuta, Erik Langerholc, and Rok Blagus. 

All the code was written by Tina Košuta and Rok Blagus. In case of questions or comments please contact rok.blagus@mf.uni-lj.si. There are no restricutions for using and/or adapting this code and or using/sharing the real data used to illustrate the methods. NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is".

This folder contains the following data and files that can be used to reproduce all analyses and figures/tables of the manuscript/Web Appendix.

It contains five subfolders:

./extensions/: [folder for reproducing the results reported in Web Appendix I - extension to non-linear mixed effects models, and mixed models with correlated random effects]
script_extensions.R: running this scripts gives all results reported in Web Appendix I. Running this script requires the following R packages: mvtnorm, nlme, pedigreemm, glmmTMB, and spaMM.

./RealData/: [folder for reproducing the results reported in Section 4.2. - Real data illustration and Web Appendix H - Additional Material for Real Data Illustration]
script_real.R: running this script performs the analysis for the 3 real data examples: the Ellenberg data (section 4.2 and Web Appendix H.I), data for assessing mathematics learning (Web Appendix H.II), and data on the number of ticks (Web Appendix H.III). The script also produces Figure 3, Figure 4, WebFigure 13, WebFigure 14, Web Table 3, Web Table 4, and Web Table 5. Running the script requires the following R packages: glmmTMB, lme4, xtable, nlme, bme, car, merDeriv, haven, tidyverse, blmeco, ggplot2, and gridExtra. The script also uses the file MathEdataset.csv, which contains the data for assessing mathematics learning 
MathEdataset.csv: the data for assessing mathematics learning used to illustrate the binomial mixed effects model.

./simulatedData/: [folder for reproducing the results reported in Section 4.1.- Simulation study and Web Appendix G - Additional Simulation Results]
./bin/: [folder containing the scripts: script_i.R, i in 1,...,40. Each script_i.R performs 1000 simulation runs for a particular simulation scenario for the binomial mixed model (e.g., the simulation settings in script_1.R are specified as: num_cluster =25, num_subj=10, sd_int=sqrt(0.01), multiplier_int=0.5, and rho=0.5), and saves the result in .txt files in the folder ./results/. The results for each method for each scenario are saved as a separate txt file. The used values of tau as determined by the proposed procedure are saved as seperate txt files for each simulation scenario. Each script uses ../source_fun.R, which actually performs all the calculations.]
./lmm/: [folder containing the scripts: script_i.R, i in 1,...,40. Each script_i.R performs 1000 simulation runs for a particular simulation scenario for the linear mixed model (e.g., the simulation settings in script_1.R are specified as: num_cluster =25, num_subj=10, sd_int=sqrt(0.01), multiplier_int=0.5, and rho=0.5), and saves the result in .txt files in the folder ./results/. The results for each method for each scenario are saved as a separate txt file. The used values of tau as determined by the proposed procedure are saved as seperate txt files for each simulation scenario.  Each script uses ../source_fun.R, which actually performs all the calculations.]
./pois/: [folder containing the scripts: script_i.R, i in 1,...,40. Each script_i.R performs 1000 simulation runs for a particular simulation scenario for the Poisson mixed model (e.g., the simulation settings in script_1.R are specified as: num_cluster =25, num_subj=10, sd_int=sqrt(0.01), multiplier_int=0.5, and rho=0.5), and saves the result in .txt files in the folder ./results/. The results for each method for each scenario are saved as a separate txt file. The used values of tau as determined by the proposed procedure are saved as seperate txt files for each simulation scenario.  Each script uses ../source_fun.R, which actually performs all the calculations.]
make_final_figs.R: this script uses the Rdata objects  res_bin.Rdata, res_lin.Rdata, and res_pois.Rdata created by read_res_bin.R, read_res_lmm.R, and read_res_pois.R, respectively and produces Figure 3 and Figure 4. It requires the ggplot2 R package.   
make_scripts.R: running this scipt produces script_i.R, i in 1,...,40 in folder ./bin/, ./lmm/, and ./pois/
read_res_bin.R: running this script reads all results (the txt files resulting from runing script_i.R, i in 1,...,40 in folder ./bin/), produces the figures reported in Web Appendix G.II, and saves the results for the binomial mixed model as res_bin.Rdata, used by  make_final_figs.R. Running this script requires the following R packages: ggplot2 and gridExtra.
read_res_lmm.R: running this script reads all results (the txt files resulting from runing script_i.R, i in 1,...,40 in folder ./lmm/), produces the figures reported in Web Appendix G.I, and saves the results for the linear mixed model as res_lin.Rdata, used by  make_final_figs.R.  Running this script requires the following R packages: ggplot2 and gridExtra.
read_res_pois.R running this script reads all results (the txt files resulting from runing script_i.R, i in 1,...,40 in folder ./pois/), produces the figures reported in Web Appendix G.III, and saves the results for the Poisson mixed model as res_pois.Rdata, used by  make_final_figs.R.
source_fun.R: this script is used by script_i.R to perform 1000 simulation runs for the specified simulation scenario: this is the core function to reproduce the simulation study. Here, the data according to the specified scenario are created and different methods are used to estimate the model's parameters. The results for each simulation scenarion and method are saved as txt files, where rows correspond to different simulation runs, and columns to the estimated parameters. Running this script requires the following R packages: glmmTMB, mvtnorm, and blme. Changing the following line of the code (L 971): alpha=0.05, to alpha=0.01 and alpha=0.1, produces the results as discussed in Section 5. We performed this part of the analysis in the cluster of CentOS based containers, session info:
R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] blme_1.0-6     lme4_1.1-34    Matrix_1.6-2   mvtnorm_1.0-11 glmmTMB_1.1.11

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.11         lattice_0.20-38     reformulas_0.4.0
 [4] rbibutils_2.3       MASS_7.3-51.4       grid_3.6.0
 [7] nlme_3.1-163        Rdpack_2.6.4        minqa_1.2.5
[10] nloptr_1.2.2.1      boot_1.3-28.1       splines_3.6.0
[13] tools_3.6.0         TMB_1.9.6           numDeriv_2016.8-1.1
[16] compiler_3.6.0      mgcv_1.9-0

