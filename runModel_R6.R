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
  #startValues$reproStat[5] <- F
lociNames <- unique(sub("[.].*$","",names(startValues)[-c(1:11)]))
genoCols = 12:ncol(startValues); startValues$age <- as.numeric(startValues$age); 

# Immigrant population
immPop <- read.csv('./Data/genotypes/immigPop.csv')
immRate <- (0/12) / 12
immMaleProb <- 1

# Demographics
surv <- read.csv('./Data/survival//survivalMonthlyREAL.csv')
ageTrans <- read.csv('./Data/stageTrans/stageTrans.csv')
probBreed <- read.csv('./Data/reproduction/probBreed_monthly.csv')

  # Choose only one of the following litterProbs
litterProbs <- read.csv('./Data/reproduction/litterProbREAL.csv')
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
iter = 1000
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
set.seed(1000)
sim1 <- simClass$new()
sim1$startParSim(numCores = numCores, iter = iter, years = years, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
              surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
              Kf = Kf, Km = Km, senesc = senesc, minMaleReproAge = minMaleReproAge, maxN_ReproMale = maxN_ReproMale,
              immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
              genOutput = genOutput, savePopulations = savePopulations, verbose = verbose)
save.image('WHATEVER.Rdata')

# Display summary statistics
sim1$summary()
sim1$plot(fieldStat=c('lambda', 'extinctTime', 'Na', 'Ho'))
#sim1$plot(fieldStat=c('pop.size', 'lambda', 'extinctTime', 'PropPoly', 'Ne', 'Na', 'Ho', 'He', 'IR', 'Fis'))

# Plot projections
matplot2(as.matrix(sim1$pop.size$All$TotalN[1:4,]))

# Pull genetic values for a given year
    yrs = c(0, 25) #, 50)
    genoMetric = c("Na",'He', 'Ho', 'PropPoly')
    sim1$pullGenoSummary(yrs, genoMetric)
    
# Pull genetic values including within population SEs
    yrs = c(0, 25)
    genoMetric = c("Na",'He', 'Ho', 'Fis', 'PropPoly')
    pullGenoWithinPop(sim1, yrs, genoMetric)

    
# Pull immigrant pop data
    # For newer sim objects
    sim1$immigrants()
    sim1$immigrants()$summary
    sim1$immigrants()$byPop
      #or
    sim1$summary()
    
# Bootstrap extinction probability
    library(boot)
    
    # bootstrapping with 1000 replications
    extinctions <- as.numeric(sim1$extinct)   #change sim1 to the name of the sim object
    extProb <- function(d, i) {
      nd <- d[i]
      return(mean(nd))
    }
    out <- boot(data=extinctions, statistic=extProb, R=50000, stype='i')
    
    # view results
    out
    plot(out)
    
    # get 95% confidence interval
    boot.ci(out, conf=0.95, type=c("basic", 'perc'))

# Bootstrap extinction probability, HPDI
    library(coda)
    
    # HPDI bootstrapping function
    bootHPD <- function(dStat, reps = 10000, prob = 0.95) {
      mSamp <- c()
      for (r in 1:reps) {
        samp <- sample(extinctions, replace=T)
        mSamp <- c(mSamp, mean(samp))
      }
      return(HPDinterval(as.mcmc(mSamp), prob = prob, na.rm = T))
    }
  
    # bootstrapping with 10000 replications
    extinctions <- as.numeric(sim1$extinct)   #change sim1 to the name of the sim object
    bootHPD(extinctions, reps=50000, prob = 0.95)
    
    
    
    
    
    
    
    
    
    
        