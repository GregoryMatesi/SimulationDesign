#install.packages("reticulate)
#library("reticulate")
source_python("/home/jovyan/work/HA_script.py")


B <- read.delim("/home/jovyan/mixtures/example_sims/Afr_CEU_10000tot_9900Afr_sims_and_reference.txt")
A <- cbind(B$afr_MAF, B$CEU_MAF)
af <- cbind(B$sample_af)
guess <- rbind(.5, .5)
HA(A, af, guess)