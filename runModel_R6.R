###########################################################################
#  Individually based model
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, John Benson

#  Publications:
# Benson, J.F., Mahoney, P.J., Sikich, J.A., Serieys, L.E., Pollinger, J.P., 
# Ernest, H.B. and Riley, S.P., 2016. Interactions between demography, genetics, 
# and landscape connectivity increase extinction probability for a small 
# population of large carnivores in a major metropolitan area. Proc. R. Soc. B, 
# 283(1837), p.20160957.

# Benson John F., Peter J. Mahoney, T. Winston Vickers, Jeff A. Sikich, 
# Paul Beier, Seth P.D. Riley, Holly B. Ernest, Walter M. Boyce.  2019.  
# Extinction vortex dynamics of top predators isolated by urbanization.  
# Ecological Applications

###########################################################################

#Load required packages
source('classes_IBM_R6.R')

##################################
## Starting genetic assignments ##
##################################
  # Starting values for SMM
startValues <- read.csv('./Data/genotypes/startValues_SMM.csv', stringsAsFactors=F)
genoCols = 12:ncol(startValues); startValues$age <- as.numeric(startValues$age);
names(startValues)[genoCols] <- gsub('[.]', '_', names(startValues)[genoCols])

  # Starting values for the SA
#startValues <- read.csv('./Data/genotypes/startValues_SA.csv', stringsAsFactors = F)
#genoCols = 15:ncol(startValues); startValues$age <- as.numeric(startValues$age);

  # Pulls unique loci names
lociNames <- unique(sub("[_].*$","",names(startValues)[genoCols]))

##################
## Demographics ##
##################
  # Survival (Choose one)
#surv <- read.csv('./Data/survival//survivalMonthly_SMM.csv')  # Santa Monicas
surv <- read.csv('./Data/survival/survivalMonthly_SA.csv')     # Santa Anas

  # Stage class transition probabilities
ageTrans <- read.csv('./Data/stageTrans/stageTrans.csv')  

  # Reproduction (Same for both populations, except for maxN_ReproMale )
probBreed <- 1        # Female probability of breeding, monthly
litterProbs <- read.csv('./Data/reproduction/litterProb.csv')  # Multinomial litter size probabilities
litterProbs$cumProbs <- cumsum(litterProbs$prob)  # Cumulative litter size probs
probFemaleKitt <- 0.5   # Probability of female kitten/offspring
senesc <- 15            # Age at senescence
minMaleReproAge <- 36   # Minimum age for reproductive males, months
maxN_ReproMale <- 5     # Reproductive K for males; 5 for Santa Anas, 2 for Santa Monicas


######################################
## Sex specific carrying capacity K ##
######################################
## (Choose one)
Kf <- matrix(c(11, 1, 0), nrow=1) ; Km <- matrix(c(5, 1, 0), nrow=1) # For SA
Kf <- matrix(c(6, 1, 0), nrow=1) ; Km <- matrix(c(2, 1, 0), nrow=1)  # For SMM

##########################
## Immigrant population ##
##########################
## Immigrant population genotypes (Choose one)
immPop <- read.csv('./Data/genotypes/immigPop_SA.csv') # SA
#immPop <- read.csv('./Data/genotypes/immigPop_SMM.csv') # SMM

# Rename genos in immPop file
names(immPop)[3:ncol(immPop)] <- gsub('[.]', '_', names(immPop)[3:ncol(immPop)])  # names of genotypes
immRate <- 0     # Immigration rate, 0 = No immigration
immMaleProb <- 1 # Probability of male immigrant, should be 0-1

# Replacement probabilities
femReplaceProb <- 0.0  # Probability a immigrant/translocated individual will replace resident
maleReplaceProb <- 0.5 # Probability a immigrant/translocated individual will replace resident



####################
## Translocations ##
####################
transMets <- NULL ; transPop <- NULL  # to turn off Translocations

##########################
## Insemination Program ##
##########################
aiRate <- NULL; aiGenotypes <- NULL  # to turn off AI




###########################
## Simulation parameters ##
###########################
genOutput <- T        # Generate genetic output, T/F
savePopulations <- T  # Save simulated populations, T/F
verbose <- T          # Provide verbose output
iter = 5000           # Number of simulated populations
years = 50            # Number of simulated years per population
numCores <- 24        # Number of cores to use if startParSim is used; detectCores() - 1


###########################
## Run model in parallel ##
###########################
set.seed(1000)
sim1 <- simClass$new()
system.time(sim1$startParSim(numCores = numCores, iter = iter, years = years, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
                 surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
                 Kf = Kf, Km = Km, senesc = senesc, minMaleReproAge = minMaleReproAge, maxN_ReproMale = maxN_ReproMale,
                 immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
                 transPop = transPop,  transMets = transMets, 
                 femReplaceProb = femReplaceProb, maleReplaceProb = maleReplaceProb,
                 aiRate = aiRate, aiGenotypes = aiGenotypes,
                 genOutput = genOutput, savePopulations = savePopulations, verbose = verbose))
save(sim1, file = 'Output.dat')


#########################
## Run model in serial ##
#########################
sim1 <- simClass$new()
sim1$startSim(iter = iter, years = years, startValues = startValues, lociNames = lociNames, genoCols = genoCols, 
              surv = surv, ageTrans = ageTrans, probBreed = probBreed, litterProbs = litterProbs, probFemaleKitt = probFemaleKitt,
              Kf = Kf, Km = Km, senesc = senesc, minMaleReproAge = minMaleReproAge, maxN_ReproMale = maxN_ReproMale,
              immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
              transPop = transPop, transMets = transMets, 
              femReplaceProb = femReplaceProb, maleReplaceProb = maleReplaceProb,
              aiRate = aiRate, aiGenotypes = aiGenotypes,
              genOutput = genOutput, savePopulations = savePopulations, verbose = verbose)
save(sim1, file = 'Output.dat')


################################
## Display summary statistics ##
################################
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
# Swap immigrants out with translocated or aiOffspring for summaries of each
sim1$immigrants()
sim1$immigrants()$summary
sim1$immigrants()$byPop
#or
sim1$summary()

