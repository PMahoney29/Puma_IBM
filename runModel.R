###########################################################################
#  Individually based model
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
###########################################################################

#Load required packages
source('classes_IBM.R')

# Functions
pullGenoSummary <- function(simObj, years, genoMetric) {
  Y <- paste("Y", years, sep="")
  o <- list()
  for (y in Y) {
    oGen <- data.frame()
    for (g in genoMetric) {
      #stat <- field(simObj, g)$mean[, y]
      stat <- simObj[[g]]$mean[,y]
      oGen <- rbind(oGen,
                    cbind(GenoMetric = g, 
                          mean = mean(stat, na.rm = T), 
                          se = sd(stat, na.rm = T) / sqrt(length(stat)),
                          lHPDI95 = HPDinterval(as.mcmc(stat), prob = 0.95, na.rm = T)[1],
                          uHPDI95 = HPDinterval(as.mcmc(stat), prob = 0.95, na.rm = T)[2]))
    }
    o[[y]] <- oGen
  }
  return(o)
}
immigrants <- function (simObj, years) {
  ps <- simObj$populations 
  #o <- aggregate(list(NumImmigrants = ps$immigrant), by = list(PopID = ps$PopID), sum)
  #o$rate <- o$NumImmigrants / years
  o <- ddply(ps, .(PopID), summarize, Immigrants = sum(immigrant), 
                                      ImmRate = sum(immigrant) / years,
                                      ReproductiveImm = sum(!is.na(reproHist) & immigrant),
                                      ReproImmRate = sum(!is.na(reproHist) & immigrant) / years)
  
  mNumImm <- mean(o$Immigrants)
  ciNumImm <- HPDinterval(as.mcmc(o$Immigrants), prob=0.95, na.rm=T) 
  mImmRate <- mean(o$ImmRate)
  ciImmRate <- HPDinterval(as.mcmc(o$ImmRate), prob=0.95, na.rm=T) 
  mNumReproImm <- mean(o$ReproductiveImm)
  ciNumReproImm <- HPDinterval(as.mcmc(o$ReproductiveImm), prob=0.95, na.rm=T) 
  mImmReproRate <- mean(o$ReproImmRate)
  ciImmReproRate <- HPDinterval(as.mcmc(o$ReproImmRate), prob=0.95, na.rm=T) 
  r1 <- cbind(mean = mNumImm, lHPDI95 = ciNumImm[1], uHPDI95 = ciNumImm[2])
  r2 <- cbind(mean = mImmRate, lHPDI95 = ciImmRate[1], uHPDI95 = ciImmRate[2])
  r3 <- cbind(mean = mNumReproImm, lHPDI95 = ciNumReproImm[1], uHPDI95 = ciNumReproImm[2])
  r4 <- cbind(mean = mImmReproRate, lHPDI95 = ciImmReproRate[1], uHPDI95 = ciImmReproRate[2])
  o.summary <- rbind(r1, r2, r3, r4) 
  row.names(o.summary) <- c("Total Immigrants", 
                            paste("Immigrant Rate (per ", years, " years)", sep=""),
                            "Total Reproductive Immigrants",
                            paste("Reproductive Immigrant Rate (per ", years, " years)", sep=""))
  
  o.list <- list(summary = o.summary, byPop = o)
  return(o.list)
}

# Starting values
startValues <- read.csv('./Data/genotypes/startValues_complete.csv', stringsAsFactors=F)
lociNames <- unique(sub("[.].*$","",names(startValues)[-c(1:11)]))
genoCols = 12:ncol(startValues); startValues$age <- as.numeric(startValues$age); 

# Immigrant population
immPop <- read.csv('./Data/genotypes/immigPop.csv')
immRate <- (1/12) / 12
immMaleProb <- 1

# Demographics
#surv <- read.csv('./Data/survival//survivalMonthly.csv')
surv <- read.csv('./Data/survival//survivalMonthlyADJUSTEDFORINBREEDINGDEPRESSION.csv')
ageTrans <- read.csv('./Data/stageTrans/stageTrans.csv')
probBreed <- read.csv('./Data/reproduction/probBreed_monthly.csv')

  # Choose only one of the following litterProbs
#litterProbs <- read.csv('./Data/reproduction/litterProb.csv')
#litterProbs$cumProbs <- cumsum(litterProbs$prob)
litterProbs <- read.csv('./Data/reproduction/PantherLitterProb.csv')
litterProbs$cumProbs <- cumsum(litterProbs$prob)

probFemaleKitt <- 0.5
senesc <- 15
#minMaleReproAge <- 42  # in months
minMaleReproAge <- ageTrans[ageTrans$sex=='M' & ageTrans$socialStat=='SubAdult', 'low']

# Sex specific carrying capacity K
# K <- c(N, prob1, prob2, ...probN); probabilities must sum to 1 across rows
Kf <- matrix(c(5, 1, 0), nrow=1)
Km <- matrix(c(2, 1, 0), nrow=1)   #Km = 2
#Km <- rbind(c(1, 0.90, 0.10), 
#            c(2, 0.50, 0.50))
#Km <- matrix(c(1, 1, 0), nrow=1)  #Km = 1


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




