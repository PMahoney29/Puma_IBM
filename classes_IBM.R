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
require(ggplot2)


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

# Calculate geometric mean
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
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
    populations = 'data.frame',
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

popClass$methods(startPop = function(startValues, ID, sex, age, mother, father, socialStat, reproStat, genoCols, genOutput) {
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
  .self$updateStats(genOutput)
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
           if (aL[[l]]$kittens[[ks]]$age < 12) {
            aL[[l]]$kittens[[ks]]$liveStat <- FALSE
            aL[[l]]$kittens[[ks]]$mortMon <- .self$time
           }
           else {
             aL[[l]]$kittens[[ks]]$socialStat <- 'SubAdult'
           }
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
popClass$methods(updateStats = function(genOutput) {
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
    
    if (genOutput) {
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
      Fis_ind <- sapply(inbreeding(g_genind, N=50), mean)
      .self$Fis[, year] <- c(mean(Fis_ind), sd(Fis_ind) / sqrt(length(Fis_ind)))
    }
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
                                     Kf, Km, genOutput = TRUE, savePopulations = TRUE, verbose = TRUE) {
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
                  socialStat='socialStat', reproStat='reproStat', genoCols=genoCols, genOutput)
    
    # simulate population
    months <- years * 12
    for (m in 1:months) {
      popi$kill(surv)
      popi$incremTime()
      popi$stageAdjust(ageTrans, Km=Km, Kf=Kf)
      popi$updateBreedStat()
      popi$reproduce(l2,l3,l4,probBreed,probFemaleKitt,lociNames)
      #popi$kill(surv)
      popi$updateStats(genOutput)
      
      if (popi$extinct == TRUE) break
    }
    
    #if (savePopulations == TRUE) field("populations", append(field("populations"), list(popi)))
    if (savePopulations == TRUE) {
      field("populations", rbind(field("populations"), cbind(PopID = popi$popID, popi$tabIndsAll())))
      #field("populations", append(field("populations"), list(popi)))
    }
    
    for (stage in 1:length(field('pop.size'))) {
      .self$pop.size[[stage]] <- rbind.fill(.self$pop.size[[stage]], popi$pop.size[stage, ])
    }
   if (genOutput) {
      for (stat in 1:length(field('Na'))) {
        .self$Na[[stat]] <- rbind.fill(.self$Na[[stat]], popi$Na[stat, ])
      }
      #for (stat in 1:length(field('Ne'))) {
    #   .self$Ne[[stat]] <- rbind.fill(.self$Ne[[stat]], popi$Ne[stat, ])
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
      field('PropPoly', rbind.fill(field('PropPoly'), popi$PropPoly))
    }
    field('lambda', rbind.fill(field('lambda'), popi$lambda))
    field('extinct', c(field('extinct'), popi$extinct))
  }
})

simClass$methods(summary = function() {
  
  N.iter <- field('iterations')
  N.years <- field('years')
  
  # Lambda
  EmpLambda_means <- sort(apply(field('lambda')[, -1], 1, function(x) gm_mean(x, na.rm=T)))
  EmpLambda_mean <- mean(EmpLambda_means)
  EmpLambda_median <- median(EmpLambda_means)
  EmpLambda_se <- sd(EmpLambda_means) / sqrt(N.iter)
  EmpLambda_quant <- quantile(EmpLambda_means, probs = c(0.025, 0.975), na.rm = T)
  
  StochLogLambda_means <- apply(field('lambda')[, -1], 1, function(x) mean(log(x), na.rm=T))
  StochLogLambda_mean <- mean(StochLogLambda_means)
  StochLogLambda_median <- median(StochLogLambda_means)
  StochLogLambda_se <- sd(StochLogLambda_means) / sqrt(N.iter)
  StochLogLambda_quant <- quantile(StochLogLambda_means, probs = c(0.025, 0.975), na.rm = T)
  

  # Extinction prob
  Prob.extinct <- mean(field('extinct'))
  
  # Mean final population size
  ps <- field('pop.size')
  outSize <- data.frame()
  for (p in 1:length(ps)) {
    ps[[p]][is.na(ps[[p]])] <- 0
    stage <- names(ps)[p]
    mean <- mean(ps[[p]][, ncol(ps[[p]])])
    se <- sd(ps[[p]][, ncol(ps[[p]])]) / sqrt(N.iter)
    outSize <- rbind(outSize, cbind(stage, mean, se))
  }
  #llply(ps, function(x) c(mean(x[, ncol(x)]), sd(x[, ncol(x)]) / sqrt(N.iter)))

  if (!is.null(.self$Na$mean)) {
    # Mean final genetics
    outGen <- data.frame()
    Nai <- .self$Na$mean[, N.years + 1]
    #Nei <- .self$Ne$mean[, N.years + 1]
    PropPolyi <- .self$PropPoly[, N.years + 1]
    Hei <- .self$He$mean[, N.years + 1]
    Hoi <- .self$Ho$mean[, N.years + 1]
    IRi <- .self$IR$mean[, N.years + 1]
    Fisi <- .self$Fis$mean[, N.years + 1]

    outGen <- rbind(outGen, 
                    cbind(stat = "Na", 
                          mean = mean(Nai, na.rm = T), 
                          se = sd(Nai, na.rm = T) / sqrt(N.iter),
                          l95_rank = c(quantile(Nai, prob = c(0.025), na.rm = T)),
                          u95_rank = c(quantile(Nai, prob = c(0.975), na.rm = T))))
    outGen <- rbind(outGen, 
                    cbind(stat = "Ne", mean = NA, se = NA, l95_rank = NA, u95_rank = NA)) 
                    #mean = mean(Nei, na.rm = T), se = sd(Nei, na.rm = T) / sqrt(N.iter)))
    outGen <- rbind(outGen, 
                    cbind(stat = "PropPoly", 
                          mean = mean(PropPolyi, na.rm = T), 
                          se = sd(PropPolyi, na.rm = T) / sqrt(N.iter),
                          l95_rank = c(quantile(PropPolyi, prob = c(0.025), na.rm = T)),
                          u95_rank = c(quantile(PropPolyi, prob = c(0.975), na.rm = T))))
    outGen <- rbind(outGen, 
                    cbind(stat = "He", 
                          mean = mean(Hei, na.rm = T), 
                          se = sd(Hei, na.rm = T) / sqrt(N.iter),
                          l95_rank = c(quantile(Hei, prob = c(0.025), na.rm = T)),
                          u95_rank = c(quantile(Hei, prob = c(0.975), na.rm = T))))
    outGen <- rbind(outGen, 
                    cbind(stat = "Ho", 
                          mean = mean(Hoi, na.rm = T), 
                          se = sd(Hoi, na.rm = T) / sqrt(N.iter),
                          l95_rank = c(quantile(Hoi, prob = c(0.025), na.rm = T)),
                          u95_rank = c(quantile(Hoi, prob = c(0.975), na.rm = T))))
    outGen <- rbind(outGen, 
                    cbind(stat = "IR", 
                          mean = mean(IRi, na.rm = T), 
                          se = sd(IRi, na.rm = T) / sqrt(N.iter),
                          l95_rank = c(quantile(IRi, prob = c(0.025), na.rm = T)),
                          u95_rank = c(quantile(IRi, prob = c(0.975), na.rm = T))))
    outGen <- rbind(outGen, 
                    cbind(stat = "Fis", 
                          mean = mean(Fisi, na.rm = T), 
                          se = sd(Fisi, na.rm = T) / sqrt(N.iter),
                          l95_rank = c(quantile(Fisi, prob = c(0.025), na.rm = T)),
                          u95_rank = c(quantile(Fisi, prob = c(0.975), na.rm = T))))
    
    row.names(outGen) <- rep(NULL, nrow(outGen))
  
    out <- list(N.iter = N.iter, N.years = N.years,
              Lambda = data.frame(mean = c(EmpLambda_mean, exp(StochLogLambda_mean)), 
                                  median = c(EmpLambda_median, exp(StochLogLambda_median)),
                                  l95_rank = c(EmpLambda_quant[1], exp(StochLogLambda_quant[1])), 
                                  u95_rank = c(EmpLambda_quant[2], exp(StochLogLambda_quant[2])),
                                  l95_se = c(EmpLambda_mean - 1.96*EmpLambda_se, exp(StochLogLambda_mean - 1.96*StochLogLambda_se)), 
                                  u95_se = c(EmpLambda_mean + 1.96*EmpLambda_se, exp(StochLogLambda_mean + 1.96*StochLogLambda_se)),
                                  row.names = c("Empirical Lambda", "Stochastic Lambda")),
              ExtinctionProb = Prob.extinct,
              Pop.size = outSize,
              GeneticComposition = outGen)
  }
  else {
    out <- list(N.iter = N.iter, N.years = N.years,
                Lambda = data.frame(mean = c(EmpLambda_mean, exp(StochLogLambda_mean)), 
                                    median = c(EmpLambda_median, exp(StochLogLambda_median)),
                                    l95_rank = c(EmpLambda_quant[1], exp(StochLogLambda_quant[1])), 
                                    u95_rank = c(EmpLambda_quant[2], exp(StochLogLambda_quant[2])),
                                    l95_se = c(EmpLambda_mean - 1.96*EmpLambda_se, exp(StochLogLambda_mean - 1.96*StochLogLambda_se)), 
                                    u95_se = c(EmpLambda_mean + 1.96*EmpLambda_se, exp(StochLogLambda_mean + 1.96*StochLogLambda_se)),
                                    row.names = c("Empirical Lambda", "Stochastic Lambda")),
                ExtinctionProb = Prob.extinct,
                Pop.size = outSize,
                GeneticComposition = "No genetic output generated")    
  }
  
  out
})

simClass$methods(plot = function(fieldStat) {
  #if (is.null(fieldStat)) fieldStat <- c('pop.size', 'lambda', 'Na', 'Ne', 'PropPoly', 'He', 'Ho', 'IR', 'Fis')
  if (is.null(fieldStat)) fieldStat <- c('pop.size', 'PropPoly', 'Na', 'He', 'Ho', 'IR', 'Fis', 'lambda')
  
  for (p in 1:length(fieldStat)) {
    par(ask=TRUE)
    
    if (fieldStat[p]=='pop.size') {
      ps <- field(fieldStat[p])
      mplots <- list()
      for (i in 1:length(ps)) {
        psi_mean <- apply(ps[[i]], 2, function(x) mean(x, na.rm=T))
        psi_low <- apply(ps[[i]], 2, function(x) quantile(x, prob = 0.025, na.rm=T))
        psi_hi <- apply(ps[[i]], 2, function(x) quantile(x, prob = 0.975, na.rm=T))
        dat_psi <-  data.frame(month = 0:(ncol(ps[[i]])-1), psi_mean = psi_mean, psi_l95 = psi_low, psi_u95 = psi_hi)
        erib <- aes(ymax = psi_u95, ymin = psi_l95)
        #psi_se <- apply(ps[[i]], 2, function(x) sd(x, na.rm=T) / sqrt(nrow(ps[[i]])))
        #dat_psi <- data.frame(month = 0:(ncol(ps[[i]])-1),psi_mean, psi_se)
        #erib <- aes(ymax = psi_mean + psi_se, ymin = psi_mean - psi_se)
        assign(paste('ps_', names(ps)[i], sep=""), ggplot(dat_psi, aes(x=month, y=psi_mean)) + geom_line(size=1.05) + geom_ribbon(erib, alpha=0.5) +
                 labs(x="Month", y=paste("Population Size:", names(ps)[i])) + #ylim(c(0,20)) + 
                 theme(axis.text.x=element_text(angle=50, size=10, vjust=0.5),
                       axis.text.y=element_text(size=10),
                       axis.title.x = element_text(size=10, vjust=-0.65),
                       axis.title.y = element_text(size=10, vjust=1)) 
        )
        mplots[[i]] <- get(paste('ps_', names(ps)[[i]], sep=""))
      }
      multiplot(plotlist=mplots, cols = 2)
    }
    
    if (fieldStat[p]=='lambda') {
      # By Emp
      lp_Empmean <- apply(field('lambda')[, -1], 1, function(x) gm_mean(x, na.rm=T))
      dat_lp <- data.frame(x = lp_Empmean)
      lp1 <- ggplot(data = dat_lp, mapping = aes(x = x)) + geom_histogram(binwidth = 0.0025) + #xlim(c(0.7, 1.1)) +
        aes(y = ..density..) + labs(x="Emp Lambda", y="Density") +
        theme(axis.title.x = element_text(size = 20, vjust = -0.65),
              axis.title.y = element_text(size = 20, vjust = 1))
    
    # By Last Year
      lp_Stochmean <- exp(apply(field('lambda')[, -1], 1, function(x) mean(log(x), na.rm=T)))
      dat_lyp <- data.frame(x = lp_Stochmean)
      lp2 <- ggplot(data = dat_lyp, mapping = aes(x = x)) + geom_histogram(binwidth = 0.0025) + #xlim(c(0.7, 1.1)) +
        aes(y = ..density..) + labs(x="Stoch Lambda", y="Density") +
        theme(axis.title.x = element_text(size = 20, vjust = -0.65),
              axis.title.y = element_text(size = 20, vjust = 1))
      multiplot(lp1, lp2, cols = 1)
      #par(ask=T) 
    }
   
    if (fieldStat[p]=='PropPoly') {
      pp <- field('PropPoly')
      pp_mean <- apply(pp, 2, function(x) mean(x, na.rm=T))
      pp_low <- apply(pp, 2, function(x) quantile(x, prob = 0.025, na.rm=T))
      pp_hi <- apply(pp, 2, function(x) quantile(x, prob = 0.975, na.rm=T))
      dat_pp <-  data.frame(year = 0:(length(pp_mean)-1), pp_mean = pp_mean, pp_l95 = pp_low, pp_u95 = pp_hi)
      erib <- aes(ymax = pp_u95, ymin = pp_l95)
      #pp_se <- apply(pp, 2, function(x) sd(x, na.rm=T) / sqrt(nrow(pp)))
      #dat_pp <- data.frame(year = 0:(length(pp_mean)-1), pp_mean = pp_mean, pp_se = pp_se)
      #erib <- aes(ymax = pp_mean + pp_se, ymin = pp_mean - pp_se)
      pp1 <- ggplot(dat_pp, aes(x=year, y=pp_mean)) + geom_line(size=1.05) + geom_ribbon(erib, alpha=0.5) +
               labs(x="Year", y="Prop of Polymorphic Loci") + ylim(c(0,1)) + 
               theme(axis.text.x=element_text(angle=50, size=20, vjust=0.5),
                     axis.text.y=element_text(size=20),
                     axis.title.x = element_text(size=20, vjust=-0.65),
                     axis.title.y = element_text(size=20, vjust=1)) 
      multiplot(pp1, cols=1)
    }
    
    if (fieldStat[p]=='Na' | fieldStat[p]=='Ne' | fieldStat[p]=='He' | fieldStat[p]=='Ho' | fieldStat[p] =='IR' | fieldStat[p]=='Fis') {
      fi <- field(fieldStat[p])$mean
      fi_mean <- apply(fi, 2, function(x) mean(x, na.rm=T))
      fi_low <- apply(fi, 2, function(x) quantile(x, prob = 0.025, na.rm=T))
      fi_hi <- apply(fi, 2, function(x) quantile(x, prob = 0.975, na.rm=T))
      dat_fi <-  data.frame(year = 0:(length(fi_mean)-1), fi_mean = fi_mean, fi_l95 = fi_low, fi_u95 = fi_hi)
      erib <- aes(ymax = fi_u95, ymin = fi_l95)
      #fi_se <- apply(fi, 2, function(x) sd(x, na.rm=T) / sqrt(nrow(fi)))
      #dat_fi <- data.frame(year = 0:(length(fi_mean)-1), fi_mean = fi_mean, fi_se = fi_se)
      #erib <- aes(ymax = fi_mean + fi_se, ymin = fi_mean - fi_se)
      fi1 <- ggplot(dat_fi, aes(x=year, y=fi_mean)) + geom_line(size=1.05) + geom_ribbon(erib, alpha=0.5) +
        labs(x="Year", y=fieldStat[p]) + #ylim(c(0,1)) + 
        theme(axis.text.x=element_text(angle=50, size=20, vjust=0.5),
              axis.text.y=element_text(size=20),
              axis.title.x = element_text(size=20, vjust=-0.65),
              axis.title.y = element_text(size=20, vjust=1)) 
      multiplot(fi1, cols=1)
    }

  }
})


