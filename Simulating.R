# "snpsampgen.R" contains a function that takes a new set of 100,000 SNP and uses those to simulate a sample population of size k.

# "SLSQPmixturesR.R" contains a function that feeds this sample population matrix into the SLSQP algorithm in R.

# "HA_script.py" contains a function that feeds this sample population matrix into the SLSQP algorithm in Python.

# At the end of every iteration of this loop will be printed the estimated pi-values, number of iterations, and time for both versions of SLSQP.

# The source scripts needed to execute this code can be found at https://github.com/GregoryMatesi/SimulationDesign

# Needs work on the 
    # outputs: proportion estimate, iters, time, seed, population, true proportions, starting guess, team member.
    # Getting the seed for the runiform draw from snpsampgen.R

# Thanks to Ian for "snpsampgen.R" and "SLSQPmixturesR.R." And to Jordan for "HA_script.py."

total.c <- read.table("/home/jovyan/resources/5ancfinalsimdat.txt")
#install.packages("nloptr")    # Needs to be installed on the server. It contains our SLSQP R function.
library("nloptr")
library("reticulate")         # Reticulate is already installed on the server. It links Python with R.

source("/home/jovyan/work/snpsampgen.R")
source("/home/jovyan/work/SLSQPmixturesR.R")
source_python("/home/jovyan/work/HA_script.py")    # source_Python() from the "reticulate" package.

START <- Sys.time()
teamMember <- "your_name"              # Put your name here.
numberSims <- 5                        # Start with 10. then try 100 or 1000.
k <-5                                  # Keep at 5. This number is input for simulating the sampling population.
numAnc <- 2                            # Set by the user. This number is only for recording in the spreadsheet.

guess1 <- 1/k                           #              

output <- c()                          # Initialize an empty output matrix.

for (i in 1:numberSims){
    
    seed.unif <- 1 # save the seed used in the uniform draw     # Seed for SNP sampling and uniform draws.
    eurfrac <- runif(1, 0.01, 0.05)    # Uniform random from user chosen interval.
    afrfrac <- 1 - eurfrac             
    easfrac <- 0
    sasfrac <- 0
    namfrac <- 0
    
    Seed.A <- 1 # Save the seed used in the sample simulation.    
    A <- snpsampgen(k, "eur_MAF" , "afr_MAF", "eas_MAF", "sas_MAF", "nam_MAF", eurfrac, afrfrac, easfrac, sasfrac, namfrac)    
    # Simulating a sample pop with new SNPs.
    
    # Calling the SLSQPmixtures function from SLSQPmixturesR.R
    Python_ <- print(SLSQPmixtures(A, k))                              # stores pi-values, iters, time.
    
    # Calling the HA function from HA_script.py                      
    af <- cbind(A$AF)                                                  # Total allele frequency vector
    A <- cbind(A$eur_MAF, A$afr_MAF, A$eas_MAF, A$sas_MAF, A$nam_MAF)  # CHANGEME: cbind(A$afr_MAF, A$CEU_MAF) etc
    
    guess <- rbind(guess1, guess1, guess1, guess1, guess1)             # kx1 vector of starting guesses.
    HA_ <- print(HA(A, af, guess))                                     # Stores pi-values, iters, time.
    
    output <- rbind(output, c(HA_[1], HA_[2], HA_[3], HA_[4], HA_[5], HA_[6], HA_[7], Python_[1], Python_[2], Python_[3], Python_[4], Python_[5], Python_[6], Python_[7], seed.unif, seed.A, teamMember, eurfrac, afrfrac, easfrac, sasfrac, namfrac))
    print(i)
    print(Sys.time() - START)
}
colnames(output) <- c("R.eur", "R.afr", "R.eas", "R.sas", "R.nam",  "R.iters", "R.time","Python.eur", "Python.afr", "Python.eas", "Python.sas", "Python.nam", "Python.iters", "Python.time", "Seed.uniform", "Seed.A", "Team.memeber", "True.eur", "True.afr", "True.eas", "True.sas", "True.nam")

END <- Sys.time()
END - START
