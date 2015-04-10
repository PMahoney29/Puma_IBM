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
#library(PopGenReport)
source('classes_IBM.R')

# Test instances of the two classes
startValues <- read.csv('./Data/genotypes/startValues_complete.csv', stringsAsFactors=F)
lociNames <- unique(sub("[.].*$","",names(startValues)[-c(1:11)]))
genoCols = 12:ncol(startValues); startValues$age <- as.numeric(startValues$age); 

#pop1 <- popClass$new(popID = 'Population_1', time=0)
#pop1$startPop(startValues=startValues, ID='animID', sex='sex', age='age', mother='mother', father='father',
              #socialStat='socialStat', reproStat='reproStat', genoCols=genoCols)

#aLit <- list(mother=pop1$indsAll[[3]], kittens = list(pop1$indsAll[[12]],pop1$indsAll[[13]]), gestation = 0)
#aLit <- append(list(aLit), list(list(mother=pop1$indsAll[[6]], kittens = list(pop1$indsAll[[10]],pop1$indsAll[[11]]), gestation = 0)))
#pop1$activeLitters <- aLit
#aL <- aLit

surv <- read.csv('./Data/survival//survivalMonthly.csv')
ageTrans <- read.csv('./Data/stageTrans/stageTrans.csv')
probBreed <- read.csv('./Data/reproduction/probBreed_monthly.csv')
litterProbs <- read.csv('./Data/reproduction/litterProb.csv')
l2 <- litterProbs[1, 'prob']
l3 <- litterProbs[2, 'prob'] + l2
l4 <- litterProbs[3, 'prob'] + l3
probFemaleKitt <- 0.5

#pop1$reproduce(l2,l3,l4,probBreed,probFemaleKitt,lociNames)

#for (i in 1:25) {
#  pop1$stageAdjust(ageTrans)
#  pop1$updateBreedStat()
#  pop1$reproduce(l2,l3,l4,probBreed,probFemaleKitt,lociNames)
#  pop1$kill(surv)
#  pop1$incremTime()
#  pop1$updateStats()
#}



