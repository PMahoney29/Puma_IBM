############################################################################
#  Individually based model
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
############################################################################

#Load required packages
library(methods)
library(adegenet)
library(plyr)
library(popbio)
library(Rhh)
library(ggplot2)
#library(PopGenReport)
source('classes_IBM.R')

# Test instances of the two classes
startValues <- read.csv('./Data/genotypes/startValues_complete.csv', stringsAsFactors=F)
lociNames <- unique(sub("[.].*$","",names(startValues)[-c(1:11)]))
genoCols = 12:ncol(startValues); startValues$age <- as.numeric(startValues$age); 

surv <- read.csv('./Data/survival//survivalMonthly.csv')
ageTrans <- read.csv('./Data/stageTrans/stageTrans.csv')
probBreed <- read.csv('./Data/reproduction/probBreed_monthly.csv')
litterProbs <- read.csv('./Data/reproduction/litterProb.csv')
l2 <- litterProbs[1, 'prob']
l3 <- litterProbs[2, 'prob'] + l2
l4 <- litterProbs[3, 'prob'] + l3
probFemaleKitt <- 0.5
Kf <- 5
Km <- 2

sim1 <- simClass$new()
sim1$startSim(iter = 5, years = 10, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
              surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
              Kf = Kf, Km = Km, genOutput = F, savePopulations = T, verbose = T)

sim1$summary()
sim1$plot(fieldStat=c('pop.size', 'lambda', 'PropPoly', 'Na', 'He', 'Ho', 'IR', 'Fis'))

