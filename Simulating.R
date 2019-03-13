# "snpsampgen.R" contains a function that takes a new set of 100,000 SNP and uses those to simulate a sample population of size k.
# "SLSQPmixturesR.R" contains a function that feeds this population matrix into the SLSQP algorithm in R.
# "HA_script.py" contains a function that feeds this sample population matrix into the SLSQP algorithm in Python.
# At the end of every iteration of this look will be printed the estimated pi-values, number of iterations, and time for both versions of SLSQP.
# The source scripts needed to execute this code can be found at https://github.com/GregoryMatesi/SimulationDesign
# Currently the most time consuming step is generating a new sample population for every iteration. Each simulation takes 8-9 second.

install.packages("nloptr")    # Needs to be installed on the server. IT contains our SLSQP R function
library("nloptr")
library("reticulate")         # Reticulate is already installed on the server. It links Python with R
source("/home/jovyan/work/snpsampgen.R")
source_python("/home/jovyan/work/HA_script.py")
source("/home/jovyan/work/SLSQPmixturesR.R")

#load(file="mixtures/Mixtures.git/total_strict_fine_maf01_atleastone.RData")

numberSims <- 100    # Start with 10. then try 100 or 1000.
k = 2
population1 ="CEU_MAF"
population2 = "afr_MAF"
popfrac1 = .5 # Uniform random from user chosen interval
popfrac2 = .5 # Uniform random from user chosen interval
guess1 = .25
guess2 = .75

for (i in 1:numberSims){
    mylist <- list()
    A <- snpsampgen(k, population1 , population2, popfrac1, popfrac2)
    # Calling the SLSQPmixtures function from SLSQPmixturesR.R
    print(SLSQPmixtures(A, k))  
    #Calling the HA function from HA_script.py                      
    af <- cbind(A$AF)    # Might just be A$af
    A <- cbind(A$CEU_MAF, A$afr_MAF)         # CHANGEME: cbind(A$afr_MAF, A$CEU_MAF) etc
    guess <- rbind(guess1, guess2)
    print(HA(A, af, guess))
    #mylist <- list(Rcode, Pythoncode)
    
}

rm(snpsampgen, SLSQPmixtures, HA, R, minimize, r, guess, af, A, i, numberSims, population1, population2, guess1, guess2, k, mylist)
