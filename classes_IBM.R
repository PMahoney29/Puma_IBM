s############################################################################
#  Individually based model classes
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
############################################################################

#Load required packages
library(methods)
library(adegenet)
library(plyr)
library(popbio)

# Test instances of the two classes
startValues <- read.csv('./Data/genotypes/startValues_complete.csv', stringsAsFactors=F)
lociNames <- unique(sub("[.].*$","",names(startValues)[-c(1:11)]))

genoCols = 12:ncol(startValues); startValues$age <- as.numeric(startValues$age); 
pop1 <- popClass$new(popID = 'Population_1', time=0)
pop1$startPop(startValues=startValues, ID='animID', sex='sex', age='age', mother='mother', father='father',
              socialStat='socialStat', reproStat='reproStat', genoCols=genoCols)

aLit <- list(mother=pop1$indsAll[[3]], kittens = list(pop1$indsAll[[12]],pop1$indsAll[[13]]), gestation = 0)
aLit <- append(list(aLit), list(list(mother=pop1$indsAll[[6]], kittens = list(pop1$indsAll[[10]],pop1$indsAll[[11]]), gestation = 0)))
pop1$activeLitters <- aLit
aL <- aLit


surv <- read.csv('./Data/survival//survivalMonthly.csv')
ageTrans <- read.csv('./Data/stageTrans/stageTrans.csv')
probBreed <- read.csv('./Data/probBreed/probBreed_monthly.csv')


###########################
##   general functions   ##
###########################

# Derive environmentally stochastic survival rates for each time step
newSurv <- function(surv) {
  ##############################################
  ## ADD CHECKS FOR NAMING OF LEVELS AND COLUMNS
  ##############################################
  out <- surv
  out[, "s"] <- unlist(lapply(1:nrow(surv), function(x) betaval(surv[x, "s"], surv[x, "se"])))
  out[, -which(names(out)=="se")]
}


#####################
## IBM classes
## popClass <- population classes that serves as a container for indClass (inividuals)
## indClass <- individual class, so far does not 'contain' popClass
#####################
popClass <- setRefClass(
  Class = 'popClass',
  fields = list(
    popID = 'character',
    indsAll = 'list',
    indsAlive = 'list',
    activeLitters = 'list',
    time = 'numeric',
    pop.size = 'numeric',
    lambda = 'numeric'
  ))

indClass <- setRefClass(
  Class = 'indClass',
  fields = list(
    animID = 'character',
    sex = 'character',
    age = 'numeric',
    mother = 'character',
    father = 'character',
    socialStat = 'character',
    reproStat = 'logical',
    reproHist = 'character',   ## another way of handling??
    liveStat = 'logical',
    birthMon = 'numeric',
    mortMon = 'numeric',
    genotype = 'data.frame'))


##########################
##   indClass methods   ##
##########################

  # Print method for individual data
indClass$methods(tab = function() {
  fields <- names(.refClassDef@fieldClasses)
  out <- data.frame()
  for (fi in fields) {
    #####  Will need to fix depending on what I decide to do with genind() object classes
    if (class(field(fi))=='data.frame') {
      #g <- c(field(fi), sep = " ")
      #out[1,fi] <- do.call(paste, g)
      out <- cbind(out, field(fi))
    }
    else {out[1,fi] <- field(fi)}
  }
  out
 })

  # Add method for including individual data to popClass object
indClass$methods(addToPop = function(popName) {
  if(class(popName)[1] != "popClass") 
    stop("Population object must be of class : 'popClass'")

  ## need to test for unique names
  
  # extending the list of individuals
  popName$indsAll <- append(popName$indsAll, .self)
})

  # Add method for female reproduction.  Number of kittens needs to be generated in advance...
indClass$methods(femBreed = function(male, numKittens, probFemaleKitt, lociNames, population) {
  if(field('sex') != "F")
    stop("Input mother is not Female")
  if(field('liveStat') != TRUE)
    stop("Input mother is not alive and cannot reproduce!")
  if(field('reproStat') != TRUE)
    stop("Input mother is not reproductive")
  
  if(male$field('sex') != "M")
    stop("Input father is not male")
  if(male$field('liveStat') != TRUE)
    stop("Input father is not alive and cannot reproduce!")
  if(male$field('reproStat') != TRUE)
    stop("Input father is not reproductive")

  pop <- length(population$indsAll)
  
  # Generate new individuals
    # Determine IDs
  idKitt <- paste('sid', seq(pop + 1, pop + numKittens), sep = "")
  
    # Determine sex
  sexKitt <- runif(numKittens, min = 0, max = 1) <= probFemaleKitt
  sexKitt <- ifelse(sexKitt==TRUE, "F", "M")
  
    # Determine genotype..completely random following Mendelian principles
  gts <- rbind(field('genotype'), male$genotype)
  genoKitt <- matrix(NA, ncol=ncol(gts), nrow=numKittens)
  for (l in 1:length(lociNames)) {
    cols <- grep(lociNames[l], names(gts))
    genoKitt[, cols] <- apply(gts[, cols], 1, function (x) sample(x, size = numKittens, replace = TRUE))
  }
  genoKitt <- as.data.frame(genoKitt)
  names(genoKitt) <- names(gts)
  
    # Determine birth month
  bm <- population$time
  
    # Loop new individuals
  for (k in 1:numKittens) {
    # gen Individual
    ind <- indClass$new(animID=idKitt[k], sex=sexKitt[k], age=0, mother=field("animID"), father=male$field("animID"), socialStat="Kitten", 
                        reproStat=FALSE, reproHist=as.character(NA), liveStat=TRUE, birthMon=bm, mortMon=as.numeric(NA), genotype=genoKitt[k,])
    
    # add to population
    ind$addToPop(population)  
  }
  
  # update individuals alive
  population$pullAlive()
  
  # Update reproHist for mother and father
  field("reproHist", paste(field("reproHist"), numKittens, sep=":"))
  male$field("reproHist", paste(male$field("reproHist"), numKittens, sep=":"))
})



##########################
##   popClass methods   ##
##########################

popClass$methods(startPop = function(startValues, ID, sex, age, mother, father, socialStat, reproStat, genoCols) {
  sv <- startValues
  for (r in 1:nrow(sv)) {
    ind <- indClass$new(animID=sv[r,ID], sex=sv[r,sex], age=sv[r,age], mother=sv[r,mother], father=sv[r,father], socialStat=sv[r,socialStat], 
                        reproStat=sv[r,reproStat], reproHist=as.character(NA), liveStat=TRUE, birthMon=as.numeric(NA), mortMon=as.numeric(NA), 
                        genotype=sv[r,genoCols])
    ind$addToPop(.self)
  }
  .self$pullAlive()
})

  # View individual data (tabulated ~ data.frame)
popClass$methods(tabIndsAll = function() {
  dat <- field('indsAll')
  out <- c()
  for (r in 1:length(dat)) {
    out = rbind(out, dat[[r]]$tab())
  }
  out
})
popClass$methods(tabAlive = function() {
  dat <- field('indsAlive')
  if (length(dat)==0) stop('No individuals are listed as alive.')
  out <- c()
  for (r in 1:length(dat)) {
    out = rbind(out, dat[[r]]$tab())
  }
  out
})

  # Pull the live individuals and store in popClass$indsAlive
popClass$methods(pullAlive = function() {
  #alive <- llply(pop1$indsAll, function(x) if (x$liveStat==TRUE) x)
  alive <- llply(field('indsAll'), function(x) if (x$liveStat==TRUE) x)
  alive <- alive[!sapply(alive, is.null)]  
  .self$indsAlive <- alive
})

  # Update population count
popClass$methods(updateCount = function() {
  if (field('time') == 0) field('pop.size', length(field('indsAlive')))
  else {field('pop.size', c(field('pop.size'), length(field('indsAlive'))))}
})

  # Update time
popClass$methods(incremTime = function() {
  field('time', field('time') + 1)
  
  # age individuals
  alive <- field("indsAlive")
  for (i in 1:length(alive)) {
    #alive[[i]]$age <- sum(alive[[i]]$age, 1, na.rm=T)
    alive[[i]]$age <- alive[[i]]$age + 1
  }  
})

  # Update lambda

  # Check breeding status
popClass$methods(updateBreedStat = function() {
  aL <- field('activeLitters')
  
  # Adjust litters based on mortality
  for (l in length(aL):1) {
    if (aL[[l]]$mother$liveStat == FALSE) {
      # kill any remaining kittens
      invisible(llply(aL[[l]]$kittens, function(x) if (x$socialStat == "Kitten") x$liveStat <- FALSE))
      
      # remove litters
      aL[[l]] <- NULL
    }
    
    else {
      # change female reproductive status after gestation
      if (aL[[l]]$gestation >= 3) {
        aL[[l]]$mother$reproStat <- TRUE  
        aL[[l]] <- NULL
      }
      
      else {
        if (length(aL[[l]]$kittens) > 0) {
          for (k in length(aL[[l]]$kittens):1) {
            # Pulling dead kittens or dispersed kittens (SubAdults)
            if (aL[[l]]$kittens[[k]]$liveStat == FALSE | aL[[l]]$kittens[[k]]$socialStat != "Kitten") aL[[l]]$kittens[[k]] <- NULL
          }
          
          #increment gestation
          if (length(aL[[l]]$kittens) == 0 & aL[[l]]$gestation < 3) aL[[l]]$gestation <- sum(aL[[l]]$gestation, 1, na.rm=T)
        }

        else {
        #increment gestation
        aL[[l]]$gestation <- sum(aL[[l]]$gestation, 1, na.rm=T)
        }
      }
    }
  }
  
  field('activeLitters', aL)
})

  # Assess reproduction
#popClass$methods(reproduce = function(probBreed) {
#  tPB <- betaval(probBreed$prob, probBreed$se)
#  alive <- field("indsAlive")
#  
#  # Assess survival
#  for (i in 1:length(alive)) {
#    ind <- alive[[i]]
#    si <- tSurv[tSurv$sex == ind$sex & tSurv$socialStat == ind$socialStat, 's']
#    ind$liveStat <- runif(1) <= si
#  }
#  
#  .self$pullAlive()
#  })

  # Assess survival
popClass$methods(kill = function(surv) {
  tSurv <- newSurv(surv)
  alive <- field("indsAlive")
  
  # Assess survival
  for (i in 1:length(alive)) {
    ind <- alive[[i]]
    si <- tSurv[tSurv$sex == ind$sex & tSurv$socialStat == ind$socialStat, 's']
    ind$liveStat <- runif(1) <= si
  }
  
  .self$pullAlive()
})

  # Asses stage transitions
popClass$methods(stageAdjust = function(ageTrans) {
  kitsAlive <- llply(field('indsAlive'), function(x) if (x$socialStat=='Kitten') x)
  subAdultAlive <- llply(field('indsAlive'), function(x) if (x$socialStat=='SubAdult') x)
  
  #kitsAlive <- llply(pop1$indsAlive, function(x) if (x$socialStat=='Kitten') x)
  kitsAlive <- kitsAlive[!sapply(kitsAlive, is.null)]
  #subAdultAlive <- llply(pop1$indsAlive, function(x) if (x$socialStat=='SubAdult') x)
  subAdultAlive <- subAdultAlive[!sapply(subAdultAlive, is.null)]
  
  for (k in 1:length(kitsAlive)) {
    kind <- kitsAlive[[k]]
    if (kind$sex == 'F') {
     if (kind$age > ageTrans[ageTrans$sex == 'F' & ageTrans$socialStat == 'Kitten', 'age']) kind$socialStat = "SubAdult" 
    }
    else {
      if (kind$age > ageTrans[ageTrans$sex == 'M' & ageTrans$socialStat == 'Kitten', 'age']) kind$socialStat = "SubAdult" 
    }
  }
  
  for (sa in 1:length(subAdultAlive)) {
    sind <- subAdultAlive[[sa]]
    if (sind$sex == 'F') {
      if (sind$age > ageTrans[ageTrans$sex == 'F' & ageTrans$socialStat == 'SubAdult', 'age']) sind$socialStat = "Adult" 
    }
    else {
      if (sind$age > ageTrans[ageTrans$sex == 'M' & ageTrans$socialStat == 'SubAdult', 'age']) sind$socialStat = "Adult" 
    }
  }
})










