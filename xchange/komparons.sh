dir24=/home/eco/work/saemix/versions/saemix2.4/saemix/R/
dir30=/home/eco/work/saemix/saemixextension/R

# Minimal changes

kompare $dir24/aaa_generics.R $dir30/aaa_generics.R

kompare $dir24/zzz.R $dir30/zzz.R

# Classes

kompare $dir24/SaemixData.R $dir30/SaemixData.R
kompare $dir24/SaemixModel.R $dir30/SaemixModel.R
kompare $dir24/SaemixRes.R $dir30/SaemixRes.R
kompare $dir24/SaemixObject.R $dir30/SaemixObject.R

# Computational functions
kompare $dir24/compute_LL.R $dir30/compute_LL.R
kompare $dir24/func_aux.R $dir30/func_aux.R
kompare $dir24/func_FIM.R $dir30/func_FIM.R

# Main algorithm
kompare $dir24/main_initialiseMainAlgo.R $dir30/main_initialiseMainAlgo.R
kompare $dir24/main_estep.R $dir30/main_estep.R
kompare $dir24/main_mstep.R $dir30/main_mstep.R
kompare $dir24/main.R $dir30/main.R

# Parameters, Simulations
kompare $dir24/func_distcond.R $dir30/func_distcond.R
kompare $dir24/func_estimParam.R $dir30/func_estimParam.R
kompare $dir24/func_simulations.R $dir30/func_simulations.R

# Plots
kompare $dir24/func_plots.R $dir30/func_plots.R

# Distcond with and without plots
kompare $dir30/func_distcond.R $dir30/func_distcond_noplot.R
