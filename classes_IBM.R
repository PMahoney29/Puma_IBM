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
    populations = 'list',
    pop.size = 'data.frame',
    lambda = 'data.frame',
    extinct = 'numeric',
    Na = 'data.frame',
    Ne = 'data.frame',
    PropPoly = 'data.frame',
    He = 'data.frame',
    Ho = 'data.frame',
    IR = 'data.frame',
    Fis = 'data.frame'
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
           aL[[l]]$kittens[[ks]]$mortMon <- .self$time
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
    #if (field('time') == 0) field('lambda', as.numeric(NA))
    #else {field('pop.size')[length(field('pop.size'))] / field('pop.size')[length(field('pop.size')) - 12]}
    
    ### Genetic metrics
    g <- pullGenos(iAlive)
    g_genind <- df2genind(createGenInput(g), sep="_")
    sink('aux')
    sumGind <- summary(g_genind)
    g_genpop <- genind2genpop(g_genind, pop = rep(field('popID'), nrow(g)))
    sink(NULL)
    
    # Allelic richness Na
    .self$Na[, year] <- c(mean(g_genind@loc.nall), sd(g_genind@loc.nall) / sqrt(length(g_genind@loc.nall)))
    #field('Na', cbind(field('Na'), 
    #                  matrix(c(mean(g_genind@loc.nall), sd(g_genind@loc.nall) / sqrt(length(g_genind@loc.nall))), ncol=1, nrow=2, byrow=F)))
    
    # Effective Ne
    
    # Proportion of alleles polymorphic
    .self$PropPoly[, year] <- mean(isPoly(g_genind, by=c("locus")))
    #field('PropPoly', c(field('PropPoly'), sum(isPoly(g_genind, by=c("locus"))) / length(g_genind@loc.names)))
    
    # Expected heterozygosity He and Ho
    Hexp <- sumGind$Hexp
    Hobs <- sumGind$Hobs
    .self$He[, year] <- c(mean(Hexp), sd(Hexp) / sqrt(length(Hexp)))
    .self$Ho[, year] <- c(mean(Hobs), sd(Hobs) / sqrt(length(Hobs)))
    #field('He', c(field('He'), matrix(c(mean(Hexp), sd(Hexp) / sqrt(length(Hexp))), ncol=1, nrow=2, byrow=F)))
    #field('Ho', c(field('Ho'), matrix(c(mean(Hobs), sd(Hobs) / sqrt(length(Hobs))), ncol=1, nrow=2, byrow=F)))
    
    # Individual heterozygosity IR
    .self$IR[, year] <- c(mean(ir(g[,-1])), sd(ir(g[,-1])) / sqrt(nrow(g)))
    #field('IR', cbind(field('IR'), matrix(c(mean(ir(g[,-1])), sd(ir(g[,-1])) / sqrt(nrow(g))), ncol=1, nrow=2, byrow=F)))
    
    # Inbreeding coefficient Fis
    Fis_ind <- sapply(inbreeding(g_genind, N=100), mean)
    .self$Fis[, year] <- c(mean(Fis_ind), sd(Fis_ind) / sqrt(length(Fis_ind)))
    #field('Fis', cbind(field('Fis'), matrix(c(mean(Fis_ind), sd(Fis_ind) / sqrt(length(Fis_ind))), ncol=1, nrow=2, byrow=F)))
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

simClass$methods(startSim = function(iter, years) {
  
  })






