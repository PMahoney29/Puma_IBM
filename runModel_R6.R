###########################################################################
#  Individually based model
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
###########################################################################

#Load required packages
source('classes_IBM_R6.R')

# Starting values
startValues <- read.csv('./Data/genotypes/startValuesFINAL.csv', stringsAsFactors=F)
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
litterProbs <- read.csv('./Data/reproduction/litterProbNEW.csv')
litterProbs$cumProbs <- cumsum(litterProbs$prob)

probFemaleKitt <- 0.5
senesc <- 15
minMaleReproAge <- 36  # in months
maxN_ReproMale <- 1

# Sex specific carrying capacity K
Kf <- matrix(c(6, 1, 0), nrow=1)
Km <- matrix(c(2, 1, 0), nrow=1)   #Km = 2


# Simulation parameters
genOutput <- T
savePopulations <- T
verbose <- T
iter = 4
years = 50
numCores <- detectCores() - 1

# Run model in serial
sim1 <- simClass$new()
sim1$startSim(iter = iter, years = years, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
              surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
              Kf = Kf, Km = Km, senesc = senesc, minMaleReproAge = minMaleReproAge, maxN_ReproMale = maxN_ReproMale,
              immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
              genOutput = genOutput, savePopulations = savePopulations, verbose = verbose)

# Run model in parallel
sim1 <- simClass$new()
sim1$startParSim(numCores = numCores, iter = iter, years = years, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
              surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
              Kf = Kf, Km = Km, senesc = senesc, minMaleReproAge = minMaleReproAge, maxN_ReproMale = maxN_ReproMale,
              immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
              genOutput = genOutput, savePopulations = savePopulations, verbose = verbose)

# Display summary statistics
sim1$summary()
sim1$plot(fieldStat=c('lambda', 'extinctTime', 'Na'))
#sim1$plot(fieldStat=c('pop.size', 'lambda', 'extinctTime', 'PropPoly', 'Ne', 'Na', 'Ho', 'He', 'IR', 'Fis'))

# Plot projections
matplot2(as.matrix(sim1$pop.size$All$TotalN[1:4,]))

# Pull genetic values for a given year
    yrs = c(0, 25) #, 50)
    genoMetric = c("Na",'He', 'Ho')
    sim1$pullGenoSummary(yrs, genoMetric)

    
# Pull immigrant pop data
    # For newer sim objects
    sim1$immigrants()
    sim1$immigrants()$summary
    sim1$immigrants()$byPop
      #or
    sim1$summary()