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

# Starting values
startValues <- read.csv('./Data/genotypes/startValues_complete.csv', stringsAsFactors=F)
lociNames <- unique(sub("[.].*$","",names(startValues)[-c(1:11)]))
genoCols = 12:ncol(startValues); startValues$age <- as.numeric(startValues$age); 

# Immigrant population
immPop <- read.csv('./Data/genotypes/immigPop.csv')
immRate <- (1/12) / 12
immMaleProb <- 1

# Demographics
surv <- read.csv('./Data/survival//survivalMonthly.csv')
ageTrans <- read.csv('./Data/stageTrans/stageTrans.csv')
probBreed <- read.csv('./Data/reproduction/probBreed_monthly.csv')

  # Choose only one of the following litterProbs
litterProbs <- read.csv('./Data/reproduction/litterProb.csv')
litterProbs$cumProbs <- cumsum(litterProbs$prob)

litterProbs <- read.csv('./Data/reproduction/PantherLitterProb.csv')
litterProbs$cumProbs <- cumsum(litterProbs$prob)

probFemaleKitt <- 0.5
senesc <- 15

# Sex specific carrying capacity K
# K <- c(N, prob1, prob2, ...probN); probabilities must sum to 1 across rows
Kf <- matrix(c(5, 1, 0), nrow=1)
Km <- rbind(c(1, 0.90, 0.10), 
            c(2, 0.50, 0.50))
#Km <- matrix(c(1, 1, 0), nrow=1)
#Km <- matrix(c(2, 1, 0), nrow=1)


# Model parameters
genOutput <- T
savePopulations <- T
verbose <- T
iter = 10
years = 25
numCores <- detectCores() - 1

# Run model in serial
sim1 <- simClass$new()
sim1$startSim(iter = iter, years = years, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
              surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
              Kf = Kf, Km = Km, senesc = senesc, 
              immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
              genOutput = genOutput, savePopulations = savePopulations, verbose = verbose)

# Run model in parallel
sim1 <- simClass$new()
sim1$startParSim(numCores = numCores, iter = iter, years = years, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
              surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
              Kf = Kf, Km = Km, senesc = senesc, 
              immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
              genOutput = genOutput, savePopulations = savePopulations, verbose = verbose)

# Display summary statistics
sim1$summary()
sim1$plot(fieldStat=c('pop.size', 'lambda', 'PropPoly', 'Na', 'Ho', 'He', 'IR', 'Fis'))

# Plot projections
matplot2(as.matrix(sim1$pop.size[[4]][1:10,]))
