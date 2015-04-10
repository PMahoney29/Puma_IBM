############################################################################
#  Individually based model classes
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
############################################################################

#Load required packages
require(methods)
require(adegenet)
require(plyr)
require(popbio)
require(Rhh)


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

# Generate litter size for COUGARS...needs to be adjusted for other species
littSize <- function(l2, l3, l4) {
  indProb <- runif(1)
  
  if (indProb <= l2) o <- 2
  if (indProb > l2 & indProb <= l3) o <- 3
  if (indProb > l3) o <- 4
  
  o
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
    lambda = 'numeric',
    extinct = 'logical'
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
  kittens <- list()
  for (k in 1:numKittens) {
    # gen Individual
    ind <- indClass$new(animID=idKitt[k], sex=sexKitt[k], age=0, mother=field("animID"), father=male$field("animID"), socialStat="Kitten", 
                        reproStat=FALSE, reproHist=as.character(NA), liveStat=TRUE, birthMon=bm, mortMon=as.numeric(NA), genotype=genoKitt[k,])

    # add to population
    ind$addToPop(population) 
    kittens <- append(kittens, ind)
  }
  
  # update individuals alive
  population$pullAlive()
  
  # update activeLitters
  population$activeLitters <- append(population$activeLitters, 
                                    list(list(mother=.self, kittens = kittens, gestation = 0)))
                                    
  # Update reproHist and reproStat for mother and father
  field("reproStat", FALSE)
  field("reproHist", paste(field("reproHist"), numKittens, sep=":"))
  male$field("reproHist", paste(male$field("reproHist"), numKittens, sep=":"))
  #male$field("reproStat", FALSE)
})



##########################
##   popClass methods   ##
##########################

popClass$methods(startPop = function(startValues, ID, sex, age, mother, father, socialStat, reproStat, genoCols) {
  sv <- startValues
  field('extinct', FALSE)
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
  field('indsAlive', alive)
})

  # Update lambda

  # Update extinction

# Assess stage transitions
popClass$methods(stageAdjust = function(ageTrans) {
  kitsAlive <- llply(field('indsAlive'), function(x) if (x$socialStat=='Kitten') x)
  subAdultAlive <- llply(field('indsAlive'), function(x) if (x$socialStat=='SubAdult') x)
  
  #kitsAlive <- llply(pop1$indsAlive, function(x) if (x$socialStat=='Kitten') x)
  kitsAlive <- kitsAlive[!sapply(kitsAlive, is.null)]
  #subAdultAlive <- llply(pop1$indsAlive, function(x) if (x$socialStat=='SubAdult') x)
  subAdultAlive <- subAdultAlive[!sapply(subAdultAlive, is.null)]
  
  if (length(kitsAlive) > 0) {
    for (k in 1:length(kitsAlive)) {
      kind <- kitsAlive[[k]]
      if (kind$sex == 'F') {
        if (kind$age > ageTrans[ageTrans$sex == 'F' & ageTrans$socialStat == 'Kitten', 'age']) kind$socialStat = "SubAdult" 
      }
      else {
        if (kind$age > ageTrans[ageTrans$sex == 'M' & ageTrans$socialStat == 'Kitten', 'age']) kind$socialStat = "SubAdult" 
      }
    }
  }
  
  if (length(subAdultAlive) > 0) {
    for (sa in 1:length(subAdultAlive)) {
      sind <- subAdultAlive[[sa]]
      if (sind$sex == 'F') {
        if (sind$age > ageTrans[ageTrans$sex == 'F' & ageTrans$socialStat == 'SubAdult', 'age']) {
          sind$socialStat = "Adult" 
          sind$reproStat = TRUE
        }
      }

      else {
        if (sind$age > ageTrans[ageTrans$sex == 'M' & ageTrans$socialStat == 'SubAdult', 'age']) {
          sind$socialStat = "Adult" 
          sind$reproStat = TRUE
        }
      }
    }
  }
})

  # Check breeding status
popClass$methods(updateBreedStat = function() {
  aL <- field('activeLitters')
  
  # Adjust litters based on mortality of mother
  if (length(aL) > 0 ) {
    for (l in length(aL):1) {
      if (aL[[l]]$mother$liveStat == FALSE) {
        # kill any remaining kittens
        if (length(aL[[l]]$kittens) > 0) {
         for (ks in 1:length(aL[[l]]$kittens)) {
           aL[[l]]$kittens[[ks]]$liveStat <- FALSE
           aL[[l]]$kittens[[ks]]$mortMonth <- .self$time
         }
        }
        #invisible(llply(aL[[l]]$kittens, function(x) if (x$socialStat == "Kitten") {
        #                                                    x$liveStat <- FALSE
        #                                                    x$mortMonth <- .self$time
        #                                                    }))
      
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
  }
})

  # Assess reproduction
popClass$methods(reproduce = function(l2,l3,l4,probBreed,probFemaleKitt,lociNames) {
  # Generate monthly probability of breeding
  tPB <- betaval(probBreed$prob, probBreed$se)
  
  # Pull reproductive adults
  f_alive <- llply(field("indsAlive"), function(x) if (x$sex=="F" & x$socialStat=="Adult" & x$reproStat==TRUE) x)
  m_alive <- llply(field("indsAlive"), function(x) if (x$sex=="M" & x$socialStat=="Adult" & x$reproStat==TRUE) x)
  #f_alive <- llply(pop1$indsAlive, function(x) if (x$sex=="F" & x$socialStat=="Adult" & x$reproStat==TRUE) x)
  #m_alive <- llply(pop1$indsAlive, function(x) if (x$sex=="M" & x$socialStat=="Adult" & x$reproStat==TRUE) x)
  f_alive <- f_alive[!sapply(f_alive, is.null)]  
  m_alive <- m_alive[!sapply(m_alive, is.null)]  
  
  #if (length(f_alive) == 0 | length(m_alive) == 0) field('extinct', TRUE)
  if (length(f_alive) != 0 & length(m_alive) != 0) {
    for (f in length(f_alive):1) {
      if (runif(1) <= tPB) {
        # Select male mate
        mate <- sample(m_alive, size = 1)
      
        # Generate number of kitts
        numKitts <- littSize(l2, l3, l4)
      
        # Breed
        #f_alive[[f]]$femBreed(mate[[1]], numKitts, probFemaleKitt, lociNames, pop1)
        f_alive[[f]]$femBreed(mate[[1]], numKitts, probFemaleKitt, lociNames, .self)
      }
    }
  }
})

  # Assess survival
popClass$methods(kill = function(surv) {
  tSurv <- newSurv(surv)
  alive <- field("indsAlive")
  
  # Assess survival
  for (i in 1:length(alive)) {
    ind <- alive[[i]]
    si <- tSurv[tSurv$sex == ind$sex & tSurv$socialStat == ind$socialStat, 's']
    ind$liveStat <- runif(1) <= si
    
    #update mortMon for individual
    if (ind$liveStat == FALSE) ind$mortMon <- .self$time
  }
  
  .self$pullAlive()
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

# Update population count
popClass$methods(updateCount = function() {
  if (field('time') == 0) field('pop.size', length(field('indsAlive')))
  else {field('pop.size', c(field('pop.size'), length(field('indsAlive'))))}
})












