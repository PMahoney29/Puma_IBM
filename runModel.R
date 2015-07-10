###########################################################################
#  Individually based model
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
###########################################################################

#Load required packages
source('classes_IBM.R')

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
#surv <- read.csv('./Data/survival//survivalMonthlyADJUSTEDFORINBREEDINGDEPRESSION.csv')
ageTrans <- read.csv('./Data/stageTrans/stageTrans.csv')
probBreed <- read.csv('./Data/reproduction/probBreed_monthly.csv')

  # Choose only one of the following litterProbs
litterProbs <- read.csv('./Data/reproduction/litterProb.csv')
litterProbs$cumProbs <- cumsum(litterProbs$prob)
#litterProbs <- read.csv('./Data/reproduction/PantherLitterProb.csv')
#litterProbs$cumProbs <- cumsum(litterProbs$prob)

probFemaleKitt <- 0.5
senesc <- 15
minMaleReproAge <- 36  # in months
#minMaleReproAge <- ageTrans[ageTrans$sex=='M' & ageTrans$socialStat=='SubAdult', 'low']

# Sex specific carrying capacity K
Kf <- matrix(c(6, 1, 0), nrow=1)
Km <- matrix(c(2, 1, 0), nrow=1)   #Km = 2


# Simulation parameters
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
              Kf = Kf, Km = Km, senesc = senesc, minMaleReproAge = minMaleReproAge,
              immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
              genOutput = genOutput, savePopulations = savePopulations, verbose = verbose)

# Run model in parallel
sim1 <- simClass$new()
sim1$startParSim(numCores = numCores, iter = iter, years = years, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
              surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
              Kf = Kf, Km = Km, senesc = senesc, minMaleReproAge = minMaleReproAge,
              immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
              genOutput = genOutput, savePopulations = savePopulations, verbose = verbose)

# Display summary statistics
sim1$summary()
sim1$plot(fieldStat=c('lambda', 'extinctTime'))
#sim1$plot(fieldStat=c('pop.size', 'lambda', 'extinctTime', 'PropPoly', 'Ne', 'Na', 'Ho', 'He', 'IR', 'Fis'))

# Plot projections
matplot2(as.matrix(sim1$pop.size$All$TotalN[1:100,]))

# Pull genetic values for a given year
    yrs = c(0, 25) #, 50)
    genoMetric = c("Na",'He', 'Ho')

    # For newer sim objects
    sim1$pullGenoSummary(yrs, genoMetric)

    # For older sim objects
    pullGenoSummary(sim1, yrs, genoMetric)
    
# Pull immigrant pop data
    # For newer sim objects
    sim1$immigrants()
    sim1$immigrants()$summary
    sim1$immigrants()$byPop
      #or
    sim1$summary()
    
    # For older sim objects
    years = 50         #whatever it is set as for the original simulation (50 in this case
    imm <- immigrants(sim1, years)
    imm$summary
    imm$byPop




