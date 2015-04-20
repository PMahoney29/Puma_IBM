############################################################################
#  Individually based model
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
############################################################################

#Load required packages
#library(methods)
#library(adegenet)
#library(plyr)
#library(popbio)
#library(Rhh)
#library(ggplot2)
#library(parallel)
#library(foreach)
#library(doParallel)
source('classes_IBM.R')

# Test instances of the two classes
startValues <- read.csv('./Data/genotypes/startValues_complete.csv', stringsAsFactors=F)
lociNames <- unique(sub("[.].*$","",names(startValues)[-c(1:11)]))
genoCols = 12:ncol(startValues); startValues$age <- as.numeric(startValues$age); 

surv <- read.csv('./Data/survival//survivalMonthly.csv')
ageTrans <- read.csv('./Data/stageTrans/stageTrans.csv')
probBreed <- read.csv('./Data/reproduction/probBreed_monthly.csv')
litterProbs <- read.csv('./Data/reproduction/litterProb.csv')
litterProbs$cumProbs <- cumsum(litterProbs$prob)
probFemaleKitt <- 0.5
Kf <- 5
Km <- 2
#Km <- rbind(c(1, 0.75), c(2, 0.25))
senesc <- 15

genOutput <- T
savePopulations <- T
verbose <- T
iter = 10
years = 25
numCores <- detectCores() - 1

sim1 <- simClass$new()
sim1$startSim(iter = iter, years = years, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
              surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
              Kf = Kf, Km = Km, senesc = senesc, genOutput = genOutput, savePopulations = savePopulations, verbose = verbose)

sim1$startParSim(numCores = numCores, iter = iter, years = years, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
              surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
              Kf = Kf, Km = Km, senesc = senesc, genOutput = genOutput, savePopulations = savePopulations, verbose = verbose)


sim1$summary()
sim1$plot(fieldStat=c('pop.size', 'lambda', 'PropPoly', 'Na', 'Ho', 'He', 'IR', 'Fis'))

# Plot projections
matplot2(as.matrix(sim1$pop.size[[4]][1:10,]))
