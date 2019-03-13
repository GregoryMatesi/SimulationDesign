
install.packages("nloptr")
library("nloptr")
library("reticulate")
source("/home/jovyan/work/snpsampgen.R")
source_python("/home/jovyan/work/HA_script.py")
source("/home/jovyan/work/SLSQPmixturesR.R")

#load(file="mixtures/Mixtures.git/total_strict_fine_maf01_atleastone.RData")

numberSims <- 100    # Start with 10. then try 100 or 1000.
k = 2
population1 ="CEU_MAF"
population2 = "afr_MAF"
guess1 = .25
guess2 = .75

for (i in 1:numberSims){
    mylist <- list()
    A <- snpsampgen(2, population1 , population2, .5, .5)
    # Calling the SLSQPmixtures function from SLSQPmixturesR.R
    #print(SLSQPmixtures(A, k))  
    #Calling the HA function from HA_script.py                      
    af <- cbind(A$AF)    # Might just be A$af
    A <- cbind(A$CEU_MAF, A$afr_MAF)         # CHANGEME: cbind(A$afr_MAF, A$CEU_MAF) etc
    guess <- rbind(guess1, guess2)
    print(HA(A, af, guess))
    #mylist <- list(Rcode, Pythoncode)
    
}

rm(snpsampgen, SLSQPmixtures, HA, R, minimize, r, guess, af, A, i, numberSims, population1, population2, guess1, guess2, k, mylist)