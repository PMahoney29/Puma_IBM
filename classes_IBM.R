############################################################################
#  Individually based model classes
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
############################################################################

### NEED TO DO:
# FIX UNIFORMLY ZERO IN INBREEDING FUNCTION (method updateStats)
# Generate Ne stat
###

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

# Pull genotypes for genetic analyses
pullGenos <- function(iAlive) {
  g <- do.call(rbind, llply(iAlive, function(x) cbind(ID = x$animID, x$genotype)))
}

createGenInput <- function(gen) {
 gi <- gen[, -1]
 cols <- seq(1, ncol(gi), by=2)

 o <- c()
 for (c in cols) {
   o <- cbind(o, paste(gi[,c], gi[,c+1], sep="_"))
 }
 o <- as.data.frame(o)
 names(o) <- names(gi)[cols]

 #o <- cbind(o, ID = gen[,1])
 o
}


#####################
## IBM classes
## simClass <- contains population iterations
## popClass <- population classes that serves as a container for indClass (individuals)
## indClass <- individual class, so far does not 'contain' popClass
#####################
simClass <- setRefClass(
  Class = 'simClass',
  fields = list(
    iterations = 'numeric',
    years = 'numeric',
    populations = 'list',
    pop.size = 'list',
    lambda = 'data.frame',
    extinct = 'numeric',
    Na = 'list',
    Ne = 'list',
    PropPoly = 'data.frame',
    He = 'list',
    Ho = 'list',
    IR = 'list',
    Fis = 'list'
  ))

popClass <- setRefClass(
  Class = 'popClass',
  fields = list(
    popID = 'character',
    indsAll = 'list',
    indsAlive = 'list',
    activeLitters = 'list',
    time = 'numeric',
    pop.size = 'data.frame',
    lambda = 'data.frame',
    extinct = 'logical',
    Na = 'data.frame',
    Ne = 'data.frame',
    PropPoly = 'data.frame',
    He = 'data.frame',
    Ho = 'data.frame',
    IR = 'data.frame',
    Fis = 'data.frame'
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
  field('pop.size', data.frame(M0 = c(0,0,0,0), row.names=c("Kittens", "SubAdults", "Adults", "Total")))
  field('Na', data.frame(Y0 = c(0,0), row.names=c("mean", "SE")))
  field('Ne', data.frame(Y0 = c(0,0), row.names=c("mean", "SE")))
  field('He', data.frame(Y0 = c(0,0), row.names=c("mean", "SE")))
  field('Ho', data.frame(Y0 = c(0,0), row.names=c("mean", "SE")))
  field('IR', data.frame(Y0 = c(0,0), row.names=c("mean", "SE")))
  field('Fis', data.frame(Y0 = c(0,0), row.names=c("mean", "SE")))
  field('lambda', data.frame(Y0 = c(0)))
  field('PropPoly', data.frame(Y0 = c(0)))
  
  # Instantiate individuals
  for (r in 1:nrow(sv)) {
    ind <- indClass$new(animID=sv[r,ID], sex=sv[r,sex], age=sv[r,age], mother=sv[r,mother], father=sv[r,father], socialStat=sv[r,socialStat], 
                        reproStat=sv[r,reproStat], reproHist=as.character(NA), liveStat=TRUE, birthMon=as.numeric(NA), mortMon=as.numeric(NA), 
                        genotype=sv[r,genoCols])
    ind$addToPop(.self)
  }
  .self$pullAlive()
  
  # Identify active litters
  iAlive <- field('indsAlive')
  aLit <- list() 
  
    # pull current kittens
  kLits <- llply(iAlive, function (x) if (x$socialStat == 'Kitten') x)
  kLits <- kLits[!sapply(kLits, is.null)]  
    # determine mothers
  mo <- unique(unlist(llply(kLits, function(x) x$mother)))
  for (m in mo) {
   mother <- llply(iAlive, function(x) if (x$animID == m) x)
   mother <- mother[!sapply(mother, is.null)]
   
   kits <- llply(iAlive, function(x) if (x$mother == m) x)
   kits <- kits[!sapply(kits, is.null)]
   
   newLitter <- list(mother = mother[[1]], kittens = kits, gestation = 0)
   aLit <- append(aLit, list(newLitter))
  }
  
  .self$activeLitters <- aLit
  
  # Update stats
  .self$updateStats()
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

# Assess stage transitions
popClass$methods(stageAdjust = function(ageTrans, Km, Kf) {
  iAlive <- field('indsAlive')
  
  kitsAlive <- llply(iAlive, function(x) if (x$socialStat=='Kitten') x)
  kitsAlive <- kitsAlive[!sapply(kitsAlive, is.null)]
  
  tsubAdultFAlive <- llply(iAlive, function(x) if (x$socialStat=='SubAdult' & x$sex == 'F' & 
                                                   x$age > ageTrans[ageTrans$sex == 'F' & ageTrans$socialStat == 'SubAdult', 'age']) x)
  tsubAdultFAlive <- tsubAdultFAlive[!sapply(tsubAdultFAlive, is.null)]
  
  tsubAdultMAlive <- llply(iAlive, function(x) if (x$socialStat=='SubAdult' & x$sex == 'M' &
                                                   x$age > ageTrans[ageTrans$sex == 'M' & ageTrans$socialStat == 'SubAdult', 'age']) x)
  tsubAdultMAlive <- tsubAdultMAlive[!sapply(tsubAdultMAlive, is.null)]
  
  adultFemalesAlive <- llply(iAlive, function(x) if (x$socialStat=='Adult' & x$sex == 'F') x)
  adultMalesAlive <- llply(iAlive, function(x) if (x$socialStat=='Adult' & x$sex == 'M') x)

  adultFemalesAlive <- adultFemalesAlive[!sapply(adultFemalesAlive, is.null)]
  adultMalesAlive <- adultFemalesAlive[!sapply(adultMalesAlive, is.null)]
  
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
  
  if (length(tsubAdultFAlive) > 0) {
    allowF <- Kf - length(adultFemalesAlive)

    if (allowF == 0) {
      invisible(llply(tsubAdultFAlive, function (x) {
       x$liveStat <- FALSE
       #x$mortMon <- popi$time
       x$mortMon <- .self$time
      }))
    }
    else {
      sampF.size <- min(length(tsubAdultFAlive), allowF)
      samp <- sample(1:length(tsubAdultFAlive), size = sampF.size)
      invisible(llply(tsubAdultFAlive[samp], function(x) {
        x$socialStat = 'Adult'
        x$reproStat = TRUE
      }))
      invisible(llply(tsubAdultFAlive[-samp], function(x) {
        x$liveStat = FALSE
        #x$mortMon <- popi$time
        x$mortMon <- .self$time
      }))
    }
  }
    
  if (length(tsubAdultMAlive) > 0) {
    allowM <- Km - length(adultMalesAlive)
    
    if (allowM == 0) {
      invisible(llply(tsubAdultMAlive, function (x) {
        x$liveStat <- FALSE
        #x$mortMon <- popi$time
        x$mortMon <- .self$time
      }))
    }
    else {
      sampM.size <- min(length(tsubAdultMAlive), allowM)
      samp <- sample(1:length(tsubAdultMAlive), size = sampM.size)
      invisible(llply(tsubAdultMAlive[samp], function(x) {
        x$socialStat = 'Adult'
        x$reproStat = TRUE
      }))
      invisible(llply(tsubAdultMAlive[-samp], function(x) {
        x$liveStat = FALSE
        #x$mortMon <- popi$time
        x$mortMon <- .self$time
      }))
    }
  }
  
  .self$pullAlive()
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
           aL[[l]]$kittens[[ks]]$mortMon <- .self$time
         }
        }
      
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

# Update population stats
popClass$methods(updateStats = function() {
  # pull living individuals
  iAlive <- field('indsAlive')
  
  # update population size
  op <- matrix(0, nrow = 4, ncol=1, byrow = F)
  op[1,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='Kitten'))))
  op[2,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='SubAdult'))))
  op[3,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='Adult'))))
  op[4,1] <- sum(op[1:3,1])
  .self$pop.size[, paste("M", field('time'), sep="")] <- op
  
  # update extinction
  if (op[4, 1] <= 1) field('extinct', TRUE)
  else {
    if (sum(unlist(llply(iAlive, function(x) sum(x$sex=='M')))) < 1) field('extinct', TRUE)
    if (sum(unlist(llply(iAlive, function(x) sum(x$sex=='F')))) < 1) field('extinct', TRUE)
  }
  
  if ((field('time')/12)%%1 == 0) {
    # Assign year for field names
    year <- paste("Y", field('time')/12, sep = "")
    
   # lambda
    if (field('time') == 0) .self$lambda[, year] <- 1
    else {.self$lambda[, year] <- .self$pop.size[nrow(field('pop.size')), ncol(field('pop.size'))] / .self$pop.size[nrow(field('pop.size')), ncol(field('pop.size')) - 12]}
    
    ### Genetic metrics
    g <- pullGenos(iAlive)
    g_genind <- df2genind(createGenInput(g), sep="_")
    
    if (Sys.info()[[1]] == 'Windows') sink('NUL')
    else {sink('aux')}
    sumGind <- summary(g_genind)
    g_genpop <- genind2genpop(g_genind, pop = rep(field('popID'), nrow(g)))
    sink(NULL)
    
    # Allelic richness Na
    .self$Na[, year] <- c(mean(g_genind@loc.nall), sd(g_genind@loc.nall) / sqrt(length(g_genind@loc.nall)))
    
    # Effective Ne
    
    # Proportion of alleles polymorphic
    .self$PropPoly[, year] <- mean(isPoly(g_genind, by=c("locus")))
    
    # Expected heterozygosity He and Ho
    Hexp <- sumGind$Hexp
    Hobs <- sumGind$Hobs
    .self$He[, year] <- c(mean(Hexp), sd(Hexp) / sqrt(length(Hexp)))
    .self$Ho[, year] <- c(mean(Hobs), sd(Hobs) / sqrt(length(Hobs)))
    
    # Individual heterozygosity IR
    .self$IR[, year] <- c(mean(ir(g[,-1])), sd(ir(g[,-1])) / sqrt(nrow(g)))
    
    # Inbreeding coefficient Fis
    Fis_ind <- sapply(inbreeding(g_genind, N=100), mean)
    .self$Fis[, year] <- c(mean(Fis_ind), sd(Fis_ind) / sqrt(length(Fis_ind)))
  }
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



  ##########################
  ##   simClass methods   ##
  ##########################

simClass$methods(startSim = function(iter, years, startValues, lociNames, genoCols, 
                                     surv, ageTrans, probBreed, litterProbs, probFemaleKitt, 
                                     Kf, Km, savePopulations = TRUE, verbose = TRUE) {
  field('iterations', iter)
  field('years', years)
  
  # Set-up list structure for outputs stats
  field('pop.size', list(kittens = c(), SubAdults = c(), Adults = c(), TotalN = c()))
  field('Na', list(mean = c(), se = c()))
  field('Ne', list(mean = c(), se = c()))
  field('He', list(mean = c(), se = c()))
  field('Ho', list(mean = c(), se = c()))
  field('IR', list(mean = c(), se = c()))
  field('Fis', list(mean = c(), se = c()))

  for (i in 1:iter) {
    if (verbose == TRUE) cat(paste('Currently on simulation: ', i, "\n", sep=""))
    # new instances of popClass
    popi <- popClass$new(popID = paste('Population_', i, sep=""), time=0)
    
    # fill with starting values
    popi$startPop(startValues=startValues, ID='animID', sex='sex', age='age', mother='mother', father='father',
                  socialStat='socialStat', reproStat='reproStat', genoCols=genoCols)
    
    # simulate population
    months <- years * 12
    for (m in 1:months) {
      popi$incremTime()
      popi$stageAdjust(ageTrans, Km=Km, Kf=Kf)
      popi$updateBreedStat()
      popi$reproduce(l2,l3,l4,probBreed,probFemaleKitt,lociNames)
      popi$kill(surv)
      popi$updateStats()
      
      if (popi$extinct == TRUE) break
    }
    
    if (savePopulations == TRUE) field("populations", list(field("populations"), list(popi)))
    
    for (stage in 1:length(field('pop.size'))) {
      .self$pop.size[[stage]] <- rbind.fill(.self$pop.size[[stage]], popi$pop.size[stage, ])
    }
    for (stat in 1:length(field('Na'))) {
      .self$Na[[stat]] <- rbind.fill(.self$Na[[stat]], popi$Na[stat, ])
    }
    #for (stat in 1:length(field('Ne'))) {
    #  .self$Ne[[stat]] <- rbind.fill(.self$Ne[[stat]], popi$Ne[stat, ])
    #}
    for (stat in 1:length(field('He'))) {
      .self$He[[stat]] <- rbind.fill(.self$He[[stat]], popi$He[stat, ])
    }
    for (stat in 1:length(field('Ho'))) {
      .self$Ho[[stat]] <- rbind.fill(.self$Ho[[stat]], popi$Ho[stat, ])
    }
    for (stat in 1:length(field('IR'))) {
      .self$IR[[stat]] <- rbind.fill(.self$IR[[stat]], popi$IR[stat, ])
    }
    for (stat in 1:length(field('Fis'))) {
      .self$Fis[[stat]] <- rbind.fill(.self$Fis[[stat]], popi$Fis[stat, ])
    }
    field('lambda', rbind.fill(field('lambda'), popi$lambda))
    field('PropPoly', rbind.fill(field('PropPoly'), popi$PropPoly))
    field('extinct', c(field('extinct'), popi$extinct))
  }
})

simClass$methods(summary = function(simObj) {
  
  N.iter <- field('iterations')
  N.years <- field('years')
  
  # Lambda
  Mean.by.year <- mean(apply(field('lambda')[, -1], 1, function(x) mean(x, na.rm=T)))
  SE.by.year <- sd(apply(field('lambda')[, -1], 1, function(x) mean(x, na.rm=T))) / sqrt(N.iter)

  L.by.LY <- c()
  for (i in 1:N.iter) {
    iL <- field('pop.size')$TotalN[i, ]
    iL <- iL[!is.na(iL)] 
    L.by.LY <- c(L.by.LY, (iL[length(iL)] / iL[1]) ^ (1/(length(iL)-1)))
  }
  Mean.by.LastYear <- mean(L.by.LY)
  SE.by.LastYear <- sd(L.by.LY) / sqrt(N.iter)
  
  #Mean.by.year <- mean(apply(sim1$field('lambda')[, -1], 1, function(x) mean(x, na.rm=T)))
  #SE.by.year <- sd(apply(sim1$field('lambda')[, -1], 1, function(x) mean(x, na.rm=T))) / sqrt(N.iter)
  
  # Extinction prob
  Prob.extinct <- mean(field('extinct'))
  
  # Mean final population size
  
  # Mean final genetics
  
})

simClass$methods(plot = function(fieldStat) {
  
})

#sim1$field('pop.size', list(kittens = c(), SubAdults = c(), Adults = c(), TotalN = c()))
#sim1$field('Na', list(mean = c(), se = c()))
#sim1$field('Ne', list(mean = c(), se = c()))
#sim1$field('He', list(mean = c(), se = c()))
#sim1$field('Ho', list(mean = c(), se = c()))
#sim1$field('IR', list(mean = c(), se = c()))
#sim1$field('Fis', list(mean = c(), se = c()))

#for (stage in 1:length(sim1$pop.size)) {
#  sim1$pop.size[[stage]] <- rbind.fill(sim1$pop.size[[stage]], popi$pop.size[stage, ])
#}
#for (stat in 1:length(sim1$field('Na'))) {
#  sim1$Na[[stat]] <- rbind.fill(sim1$Na[[stat]], popi$Na[stat, ])
#}
#for (stat in 1:length(sim1$field('Ne'))) {
#  sim1$Ne[[stat]] <- rbind.fill(sim1$Ne[[stat]], popi$Ne[stat, ])
#}
#for (stat in 1:length(sim1$field('He'))) {
#  sim1$He[[stat]] <- rbind.fill(sim1$He[[stat]], popi$He[stat, ])
#}
#for (stat in 1:length(sim1$field('Ho'))) {
#  sim1$Ho[[stat]] <- rbind.fill(sim1$Ho[[stat]], popi$Ho[stat, ])
#}
#for (stat in 1:length(sim1$field('IR'))) {
#  sim1$IR[[stat]] <- rbind.fill(sim1$IR[[stat]], popi$IR[stat, ])
#}
#for (stat in 1:length(sim1$field('Fis'))) {
#  sim1$Fis[[stat]] <- rbind.fill(sim1$Fis[[stat]], popi$Fis[stat, ])
#} 
#sim1$field('lambda', rbind.fill(sim1$field('lambda'), popi$lambda))
#sim1$field('PropPoly', rbind.fill(sim1$field('PropPoly'), popi$PropPoly))
#sim1$field('extinct', c(sim1$field('extinct'), popi$extinct))



