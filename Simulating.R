# "snpsampgen.R" contains a function that takes a new set of 100,000 SNP and uses those to simulate a sample population of size k.

# "SLSQPmixturesR.R" contains a function that feeds this population matrix into the SLSQP algorithm in R.

# "HA_script.py" contains a function that feeds this sample population matrix into the SLSQP algorithm in Python.

# At the end of every iteration of this look will be printed the estimated pi-values, number of iterations, and time for both versions of SLSQP.

# The source scripts needed to execute this code can be found at https://github.com/GregoryMatesi/SimulationDesign

# Needs work on the 
    # outputs: proportion estimate, iters, time, seed, population, true proportions, starting guess, team member.
    # Getting the seed for the runiform draw from snpsampgen.R

# Thanks to Ian for "snpsampgen.R" and "SLSQPmixturesR.R." And to Jordan for "HA_script.py"


install.packages("nloptr")    # Needs to be installed on the server. IT contains our SLSQP R function
library("nloptr")
library("reticulate")         # Reticulate is already installed on the server. It links Python with R

rm(snpsampgen, SLSQPmixtures, HA, R, minimize, r, guess, af, A, i, numberSims, population1, population2, guess1, guess2, k, popfrac1, popfrac2, END, START, HA_, Python_) # remove everything except for total.c because it is a huge file.

source("/home/jovyan/work/snpsampgen.R")
source("/home/jovyan/work/SLSQPmixturesR.R")
source_python("/home/jovyan/work/HA_script.py")    # From the "reticulate" package.

#load(file="mixtures/Mixtures.git/total_strict_fine_maf01_atleastone.RData")

START <- Sys.time()
numberSims <- 100              # Start with 10. then try 100 or 1000.
k <- 2                       # Number of ancestries.
population1 <- "CEU_MAF"     # Name your ancestries
population2 <- "afr_MAF"     # 
popfrac1 <- runif(1, 0, 0.5) # Uniform random from user chosen interval.
popfrac2 <- 1 - popfrac1     # Uniform random from user chosen interval.
guess1 <- .25                # User chosen starting value.
guess2 <- 1 - guess1         # User chosen sarting value.
output <- c()                # initialize an empty output matrix.
teamMember <- "your_name"    # Put your name here.
seed <- 1                    # get this from snpsampgen.R. the same seed is used for both python and r versions of SLSQP.

for (i in 1:numberSims){
    
    A <- snpsampgen(k, population1 , population2, popfrac1, popfrac2)
    # Calling the SLSQPmixtures function from SLSQPmixturesR.R
    Python_ <- print(SLSQPmixtures(A, k))  
    #Calling the HA function from HA_script.py                      
    af <- cbind(A$AF)    # Might just be A$af
    A <- cbind(A$CEU_MAF, A$afr_MAF)         # CHANGEME: cbind(A$afr_MAF, A$CEU_MAF) etc
    guess <- rbind(guess1, guess2)
    HA_ <- print(HA(A, af, guess))
    output <- rbind(output, c(HA_[1], HA_[2], HA_[3], Python_[1], Python_[2], Python_[3], seed, population1, popfrac1, guess1, population2, popfrac2, guess2, teamMember))
    print(i)
}
END <- Sys.time()
END - START
