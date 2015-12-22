############################################################################
#  Individually based model classes
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
############################################################################

### NEED TO DO:
###
###

#Load required packages
require(methods)
require(adegenet)
require(plyr)
require(popbio)
require(Rhh)
require(ggplot2)
require(coda)
require(parallel)
require(foreach)
require(doParallel)
require(R6)
filename <- 'classes_IBM_R6.R'


###########################
##   general functions   ##
###########################

# Combine output from multiple simulation runs
combineSims <- function(simObjects) {
  s_out <- simClass$new()
  #s_out <<- simObjects[[1]]
  #s_out$populations$PopID <- paste('Run1_', s_out$populations$PopID, sep='')
  for (s in 1:length(simObjects)) {
    if (s==1) {
      s_out$Date <- simObjects[[s]]$Date
      s_out$SimTime <- simObjects[[s]]$SimTime
      s_out$iterations <- simObjects[[s]]$iterations
      s_out$years <- simObjects[[s]]$years
      
      simObjects[[s]]$populations$PopID <- paste('Run', s, '_', simObjects[[s]]$populations$PopID, sep='')
      s_out$populations <- simObjects[[s]]$populations
      
      s_out$pop.size <- simObjects[[s]]$pop.size
      s_out$lambda <- simObjects[[s]]$lambda
      s_out$extinct <- simObjects[[s]]$extinct
      s_out$extinctTime <- simObjects[[s]]$extinctTime
      s_out$Na <- simObjects[[s]]$Na
      s_out$Ne <- simObjects[[s]]$Ne
      s_out$PropPoly <- simObjects[[s]]$PropPoly
      s_out$He <- simObjects[[s]]$He
      s_out$Ho <- simObjects[[s]]$Ho
      s_out$IR <- simObjects[[s]]$IR
      s_out$Fis <- simObjects[[s]]$Fis
    }
    else {
      s_out$Date <- c(s_out$Date, simObjects[[s]]$Date)
      s_out$SimTime <- c(s_out$SimTime, simObjects[[s]]$SimTime)
      s_out$iterations <- c(s_out$iterations, simObjects[[s]]$iterations)

      simObjects[[s]]$populations$PopID <- paste('Run', s, '_', simObjects[[s]]$populations$PopID, sep='')
      s_out$populations <- rbind(s_out$populations, simObjects[[s]]$populations)
      
      if(s_out$years[1] != simObjects[[s]]$years) stop(simpleError('Projection lengths do not match!'))
      for (sex in 1:length(s_out$pop.size)) {
        for (stage in 1:length(s_out$pop.size[[sex]])) {
          s_out$pop.size[[sex]][[stage]] <- rbind(s_out$pop.size[[sex]][[stage]], 
                                                simObjects[[s]]$pop.size[[sex]][[stage]])
        }
      }
      s_out$lambda <- rbind(s_out$lambda, simObjects[[s]]$lambda)
      s_out$extinct <- c(s_out$extinct, simObjects[[s]]$extinct)
      s_out$extinctTime <- c(s_out$extinctTime, simObjects[[s]]$extinctTime)
      
      for (na in 1:2) {
        s_out$Na[[na]] <- rbind(s_out$Na[[na]], simObjects[[s]]$Na[[na]])
      }
      
      s_out$Ne <- rbind(s_out$Ne, simObjects[[s]]$Ne)
      
      s_out$PropPoly <- rbind(s_out$PropPoly, simObjects[[s]]$PropPoly) 
      
      for (he in 1:2) {
        s_out$He[[he]] <- rbind(s_out$He[[he]], simObjects[[s]]$He[[he]])
      } 
      
      for (ho in 1:2) {
        s_out$Ho[[ho]] <- rbind(s_out$Ho[[ho]], simObjects[[s]]$Ho[[ho]])
      } 
      
      for (ir in 1:2) {
        s_out$IR[[ir]] <- rbind(s_out$IR[[ir]], simObjects[[s]]$IR[[ir]])
      } 
      
      for (fis in 1:2) {
        s_out$Fis[[fis]] <- rbind(s_out$Fis[[fis]], simObjects[[s]]$Fis[[fis]])
      } 
    }

  }
  s_out$iterations <- sum(s_out$iterations)
  return(s_out)
}

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
littSize <- function(litterProbs) {
  if (signif(sum(litterProbs$prob)) != 1) 
    stop("Litter size probabilities must sum to 1")
  indProb <- runif(1)
  high <- min(which(litterProbs$cumProbs > indProb))
  o <- litterProbs[high, 'LitterSize']
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
simClass <- R6Class('simClass',
  public = list(
    Date = NA,    #Change to appropriate posix class
    SimTime = NA, #Change to appropriate posix class
    iterations = NA,
    years = NA,
    populations = data.frame(),
    pop.size = list(),
    lambda = data.frame(),
    extinct = FALSE,
    extinctTime = NA,  #Change to appropriate numeric with na class
    Na = list(),
    Ne = data.frame(),
    PropPoly = data.frame(),
    He = list(),
    Ho = list(),
    IR = list(),
    Fis = list(),
    
    startSim = function(iter, years, startValues, lociNames, genoCols, 
                        surv, ageTrans, probBreed, litterProbs, probFemaleKitt, 
                        Kf, Km, senesc, minMaleReproAge, maxN_ReproMale = maxN_ReproMale,
                        immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
                        genOutput = TRUE, savePopulations = TRUE, verbose = TRUE) {
      start <- Sys.time()
      self$Date <- start
      self$iterations <- iter
      self$years <- years
      immPop_subset <- immPop
      
      # Set-up list structure for outputs stats
      self$pop.size <- list(Females = list(kittens = c(), SubAdults = c(), Adults = c(), TotalN = c()),
                             Males = list(kittens = c(), SubAdults = c(), Adults = c(), TotalN = c()),
                             All = list(kittens = c(), SubAdults = c(), Adults = c(), TotalN = c()))
      self$Na <- list(mean = c(), se = c())
      self$He <- list(mean = c(), se = c())
      self$Ho <- list(mean = c(), se = c())
      self$IR <- list(mean = c(), se = c())
      self$Fis <- list(mean = c(), se = c())
      
      for (i in 1:iter) {
        if (verbose == TRUE) cat(paste('Currently on simulation: ', i, "\n", sep=""))
        
        # new instances of popClass
        popi <- popClass$new(popID = paste('Population_', i, sep=""), time=0)
        
        # fill with starting values
        popi$startPop(startValues=startValues, ID='animID', sex='sex', age='age', mother='mother', father='father',
                      socialStat='socialStat', reproStat='reproStat', genoCols=genoCols, genOutput=genOutput, ageTrans = ageTrans)
        
        # simulate population
        months <- years * 12
       # ent = Sys.time()
        for (m in 1:months) {
          print(paste('Step:', m/12))
          popi$kill(surv)
          popi$incremTime(senesc)
          if (runif(1) <= immRate) {
            immPop_subset <- popi$addImmigrants(iP=immPop_subset, immMaleProb, ageTrans)
            if (nrow(immPop_subset) == 0) immPop_subset <- immPop
          }
        #  start = Sys.time()
          popi$stageAdjust(ageTrans, Km=Km, Kf=Kf, minMaleReproAge=minMaleReproAge)
        #  print(paste('Computation time (Stage Adjust): ', (Sys.time() - start), sep=''))
          
        #  start = Sys.time()
          popi$updateBreedStat(ageTrans, maxN_ReproMale)
        #  print(paste('Computation time (Breed Status): ', (Sys.time() - start), sep=''))
          
        #  start = Sys.time()
          popi$reproduce(litterProbs,probBreed,probFemaleKitt,lociNames)
        #  print(paste('Computation time (Reproduce): ', (Sys.time() - start), sep=''))
          
        #  start = Sys.time()
          popi$updateStats(genOutput)
        #  print(paste('Computation time (Update Stats): ', (Sys.time() - start), sep=''))
          
          
          if (popi$extinct == TRUE) break
        #  print(paste('Total Time: ', Sys.time() - ent, sep=''))
        }
        
        if (savePopulations == TRUE) 
          self$populations <- rbind(self$populations, cbind(PopID = popi$popID, popi$tabAll()))

        for (subP in 1:length(self$pop.size)) {
          for (stage in 1:length(self$pop.size[[subP]])) {
            self$pop.size[[subP]][[stage]] <- rbind.fill(self$pop.size[[subP]][[stage]], popi$pop.size[[subP]][stage, ])
          }
        }
        
        if (genOutput) {
          for (stat in 1:length(self$Na)) {
            self$Na[[stat]] <- rbind.fill(self$Na[[stat]], popi$Na[stat, ])
          }
          for (stat in 1:length(self$He)) {
            self$He[[stat]] <- rbind.fill(self$He[[stat]], popi$He[stat, ])
          }
          for (stat in 1:length(self$Ho)) {
            self$Ho[[stat]] <- rbind.fill(self$Ho[[stat]], popi$Ho[stat, ])
          }
          for (stat in 1:length(self$IR)) {
            self$IR[[stat]] <- rbind.fill(self$IR[[stat]], popi$IR[stat, ])
          }
          for (stat in 1:length(self$Fis)) {
            self$Fis[[stat]] <- rbind.fill(self$Fis[[stat]], popi$Fis[stat, ])
          }
          self$PropPoly <- rbind.fill(self$PropPoly, popi$PropPoly)
          self$Ne <- rbind.fill(self$Ne, popi$Ne)
        }
        self$lambda <- rbind.fill(self$lambda, popi$lambda)
        self$extinct <- c(self$extinct, popi$extinct)
        if (popi$extinct)
          self$extinctTime <- c(self$extinctTime, (popi$time / 12))
      }
      
      self$SimTime <- (Sys.time() - start)
    },
    
    startParSim = function(numCores = detectCores(), iter, years, startValues, lociNames, genoCols, 
                           surv, ageTrans, probBreed, litterProbs, probFemaleKitt, 
                           Kf, Km, senesc, minMaleReproAge,maxN_ReproMale = maxN_ReproMale,
                           immPop = immPop, immRate = immRate, immMaleProb = immMaleProb,
                           genOutput = TRUE, savePopulations = TRUE, verbose = TRUE) {
      start <- Sys.time()
      self$Date <- start
      self$iterations <- iter
      self$years <- years
      immPop_subset <- immPop
      
      cSim <- function(s1, s2) {
        s1$populations <- rbind(s1$populations, s2$populations)
        for (subP in 1:length(s1$pop.size)) {
          for (stage in 1:length(s1$pop.size[[subP]])) {
            s1$pop.size[[subP]][[stage]] <- rbind.fill(s1$pop.size[[subP]][[stage]], s2$pop.size[[subP]][[stage]])
          }
        }
        
        s1$lambda <- rbind.fill(s1$lambda, s2$lambda)
        s1$extinct <- c(s1$extinct, s2$extinct)
        s1$extinctTime <- c(s1$extinctTime, s2$extinctTime)
        
        for (stat in 1:length(s1$Na)) {
          s1$Na[[stat]] <- rbind.fill(s1$Na[[stat]], s2$Na[[stat]])
        }
        
        for (stat in 1:length(s1$He)) {
          s1$He[[stat]] <- rbind.fill(s1$He[[stat]], s2$He[[stat]])
        }
        for (stat in 1:length(s1$Ho)) {
          s1$Ho[[stat]] <- rbind.fill(s1$Ho[[stat]], s2$Ho[[stat]])
        }
        for (stat in 1:length(s1$IR)) {
          s1$IR[[stat]] <- rbind.fill(s1$IR[[stat]], s2$IR[[stat]])
        }
        for (stat in 1:length(s1$Fis)) {
          s1$Fis[[stat]] <- rbind.fill(s1$Fis[[stat]], s2$Fis[[stat]])
        }
        s1$PropPoly <- rbind.fill(s1$PropPoly, s2$PropPoly)
        s1$Ne <- rbind.fill(s1$Ne, s2$Ne)
        
        return(s1)
      }
      
      packList <- c('methods','plyr', 'popbio', 'Rhh')
      
      #start cluster
      cl <- makeCluster(numCores)
      registerDoParallel(cl, cores = numCores)
      
      #operation
      o <- foreach(i = 1:iter, .combine = cSim, 
                   .packages = packList, 
                   .inorder = F, .verbose = FALSE, .export = 'filename') %dopar% {
                     source(filename)
                     
                     # new instances of popClass
                     popi <- popClass$new(popID = paste('Population_', i, sep=""), time=0)
                     
                     # fill with starting values
                     popi$startPop(startValues=startValues, ID='animID', sex='sex', age='age', mother='mother', father='father',
                                   socialStat='socialStat', reproStat='reproStat', genoCols=genoCols, genOutput=genOutput, ageTrans=ageTrans)
                     
                     # simulate population
                     months <- years * 12
                     for (m in 1:months) {
                       popi$kill(surv)
                       popi$incremTime(senesc)
                       if (runif(1) <= immRate) {
                         immPop_subset <- popi$addImmigrants(iP=immPop_subset, immMaleProb, ageTrans)
                         if (nrow(immPop_subset) == 0) immPop_subset <- immPop
                       }
                       popi$stageAdjust(ageTrans, Km=Km, Kf=Kf, minMaleReproAge=minMaleReproAge)
                       popi$updateBreedStat(ageTrans, maxN_ReproMale)
                       popi$reproduce(litterProbs,probBreed,probFemaleKitt,lociNames)
                       #popi$kill(surv)
                       popi$updateStats(genOutput)
                       if (popi$extinct == TRUE) break
                     }
                     
                     out <- list()
                     if (savePopulations == TRUE) {
                       out$populations <- cbind(PopID = popi$popID, popi$tabAll())
                     }
                     else {out$populations <- NULL}
                     
                     out$pop.size <- llply(popi$pop.size, function(x) {
                       list(Kittens = x[1,], SubAdults = x[2,], Adults = x[3,], TotalN = x[4,])
                     })
                     
                     out$lambda <- popi$lambda
                     out$extinct <- popi$extinct
                     if (out$extinct)
                       out$extinctTime <- popi$time / 12
                     else {out$extinctTime <- NA_real_}
                     
                     if (genOutput) {
                       out$Na$mean <- popi$Na[1, ]
                       out$Na$se <- popi$Na[2, ]
                       
                       out$He$mean <- popi$He[1, ]
                       out$He$se <- popi$He[2, ]              
                       
                       out$Ho$mean <- popi$Ho[1, ]
                       out$Ho$se <- popi$Ho[2, ]              
                       
                       out$IR$mean <- popi$IR[1, ]
                       out$IR$se <- popi$IR[2, ]
                       
                       out$Fis$mean <- popi$Fis[1, ]
                       out$Fis$se <- popi$Fis[2, ]
                       
                       out$PropPoly <- popi$PropPoly
                       out$Ne <- popi$Ne
                     }
                     else {
                       out$Na$mean <- NULL
                       out$Na$se <- NULL
                       
                       out$He$mean <- NULL
                       out$He$se <- NULL              
                       
                       out$Ho$mean <- NULL
                       out$Ho$se <- NULL              
                       
                       out$IR$mean <- NULL
                       out$IR$se <- NULL
                       
                       out$Fis$mean <- NULL
                       out$Fis$se <- NULL
                       
                       out$PropPoly <- NULL
                       out$Ne <- NULL
                     }
                     return(out)
                   }
      
      #stop cluster
      stopCluster(cl)
      
      # Set-up list structure for outputs stats
      self$populations <- o$populations
      self$pop.size <- o$pop.size
      self$Na <- list(mean = o$Na$mean, se = o$Na$se)
      self$Ne <- o$Ne
      self$He <- list(mean = o$He$mean, se = o$He$se)
      self$Ho <- list(mean = o$Ho$mean, se = o$Ho$se)
      self$IR <- list(mean = o$IR$mean, se = o$IR$se)
      self$Fis <- list(mean = o$Fis$mean, se = o$Fis$se)
      self$PropPoly <- o$PropPoly
      self$lambda <- o$lambda
      if (!is.null(o$extinct)) self$extinct <- o$extinct
      self$extinctTime <- o$extinctTime
      
      self$SimTime <- (Sys.time() - start)
      #print(paste('Computation time: ', (Sys.time() - start), sep=''))
    },
    
    summary = function() {
      
      N.iter <- self$iterations
      N.years <- self$years
      
      # Lambda
      EmpLambda_means <- sort(apply(self$lambda[, -1], 1, function(x) gm_mean(x, na.rm=T)))
      EmpLambda_mean <- mean(EmpLambda_means)
      EmpLambda_median <- median(EmpLambda_means)
      EmpLambda_se <- sd(EmpLambda_means) / sqrt(N.iter)
      EmpLambda_quant <- HPDinterval(as.mcmc(EmpLambda_means), prob = 0.95, na.rm = T)
      
      StochLogLambda_means <- apply(self$lambda[, -1], 1, function(x) mean(log(x), na.rm=T))
      StochLogLambda_mean <- mean(StochLogLambda_means)
      StochLogLambda_median <- median(StochLogLambda_means)
      StochLogLambda_se <- sd(StochLogLambda_means) / sqrt(N.iter)
      StochLogLambda_quant <- HPDinterval(as.mcmc(StochLogLambda_means), prob = 0.95, na.rm = T)
      
      # Extinction prob
      Prob.extinct <- mean(self$extinct)
      Extinct.time_mean <- mean(self$extinctTime, na.rm=T)
      Extinct.time_median <- median(self$extinctTime, na.rm=T)
      if (length(na.omit(self$extinctTime)) > 1) {
        Extinct.time_quant <- HPDinterval(as.mcmc(self$extinctTime), prob = 0.95, na.rm = T)
        eTime <- data.frame(mean = Extinct.time_mean,
                            median = Extinct.time_median,
                            lHPDI95 = Extinct.time_quant[1],
                            uHPDI95 = Extinct.time_quant[2])}
      else {
        eTime <- data.frame(mean = Extinct.time_mean,
                            median = Extinct.time_median,
                            lHPDI95 = 'Inestimable',
                            uHPDI95 = 'Inestimable') 
      }
      
      # Mean final population size
      if (((N.years * 12) + 1) == ncol(self$pop.size$All$TotalN)) {
        ps <- self$pop.size
        outSize <- llply(ps, function (x) {
          ldply(x, function (y) {
            y[is.na(y)] <- 0
            mean <- mean(y[, ncol(y)])
            se <- sd(y[, ncol(y)]) / sqrt(N.iter)
            HPDI95 = HPDinterval(as.mcmc(y[, ncol(y)]), prob = 0.95, na.rm = T)
            return(cbind(mean, se, lHPDI95 = HPDI95[1], uHPDI95 = HPDI95[2]))
          })
        })
      }
      else {outSize <- 'Final Population size = 0 since ExtinctionProb = 1'}
      
      # Immigrant population stats
      pps <- self$populations
      o <- ddply(pps, .(PopID), summarize, Immigrants = sum(immigrant), 
                 ImmRate = sum(immigrant) / years,
                 ReproductiveImm = sum(!is.na(reproHist) & immigrant),
                 ReproImmRate = sum(!is.na(reproHist) & immigrant) / years)
      
      if(!is.null(o)) {
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
        imm <- rbind(r1, r2, r3, r4) 
        row.names(imm) <- c("Total Immigrants", 
                            paste("Immigrant Rate (per ", years, " years)", sep=""),
                            "Total Reproductive Immigrants",
                            paste("Reproductive Immigrant Rate (per ", years, " years)", sep=""))
      }
      else {
        r1 <- 'N/A'
        r2 <- 0
        r3 <- 'N/A'
        r4 <- 'N/A'
        
        imm <- rbind(r1, r2, r3, r4) 
        row.names(imm) <- c("Total Immigrants", 
                            paste("Immigrant Rate (per ", years, " years)", sep=""),
                            "Total Reproductive Immigrants",
                            paste("Reproductive Immigrant Rate (per ", years, " years)", sep=""))
      }
      
      #if (!is.null(self$Na$mean) & (N.years + 1) == ncol(self$Na$mean)) {
      if (!is.null(self$Na$mean)) {
        if ((N.years + 1) == ncol(self$Na$mean)) {
      
          # Mean final genetics
          outGen <- data.frame()
          Nai <- self$Na$mean[, N.years + 1]
          Nei <- self$Ne[, N.years + 1]
          PropPolyi <- self$PropPoly[, N.years + 1]
          Hei <- self$He$mean[, N.years + 1]
          Hoi <- self$Ho$mean[, N.years + 1]
          IRi <- self$IR$mean[, N.years + 1]
          Fisi <- self$Fis$mean[, N.years + 1]
        
          if (length(na.omit(Nai)) > 1) {
            outGen <- rbind(outGen, 
                            cbind(stat = "Na", 
                                  mean = mean(Nai, na.rm = T), 
                                  se = sd(Nai, na.rm = T) / sqrt(N.iter),
                                  lHPDI95 = HPDinterval(as.mcmc(Nai), prob = 0.95, na.rm = T)[1],
                                  uHPDI95 = HPDinterval(as.mcmc(Nai), prob = 0.95, na.rm = T)[2]))
            outGen <- rbind(outGen, 
                            cbind(stat = "Ne", 
                                  mean = mean(Nei, na.rm = T), 
                                  se = sd(Nei, na.rm = T) / sqrt(N.iter),
                                  lHPDI95 = HPDinterval(as.mcmc(Nei), prob = 0.95, na.rm = T)[1],
                                  uHPDI95 = HPDinterval(as.mcmc(Nei), prob = 0.95, na.rm = T)[2]))
            outGen <- rbind(outGen, 
                            cbind(stat = "PropPoly", 
                                  mean = mean(PropPolyi, na.rm = T), 
                                  se = sd(PropPolyi, na.rm = T) / sqrt(N.iter),
                                  lHPDI95 = HPDinterval(as.mcmc(PropPolyi), prob = 0.95, na.rm = T)[1],
                                  uHPDI95 = HPDinterval(as.mcmc(PropPolyi), prob = 0.95, na.rm = T)[2]))
            outGen <- rbind(outGen, 
                            cbind(stat = "He", 
                                  mean = mean(Hei, na.rm = T), 
                                  se = sd(Hei, na.rm = T) / sqrt(N.iter),
                                  lHPDI95 = HPDinterval(as.mcmc(Hei), prob = 0.95, na.rm = T)[1],
                                  uHPDI95 = HPDinterval(as.mcmc(Hei), prob = 0.95, na.rm = T)[2]))
            outGen <- rbind(outGen, 
                            cbind(stat = "Ho", 
                                  mean = mean(Hoi, na.rm = T), 
                                  se = sd(Hoi, na.rm = T) / sqrt(N.iter),
                                  lHPDI95 = HPDinterval(as.mcmc(Hoi), prob = 0.95, na.rm = T)[1],
                                  uHPDI95 = HPDinterval(as.mcmc(Hoi), prob = 0.95, na.rm = T)[2]))
            outGen <- rbind(outGen, 
                            cbind(stat = "IR", 
                                  mean = mean(IRi, na.rm = T), 
                                  se = sd(IRi, na.rm = T) / sqrt(N.iter),
                                  lHPDI95 = HPDinterval(as.mcmc(IRi), prob = 0.95, na.rm = T)[1],
                                  uHPDI95 = HPDinterval(as.mcmc(IRi), prob = 0.95, na.rm = T)[2]))
            outGen <- rbind(outGen, 
                            cbind(stat = "Fis", 
                                  mean = mean(Fisi, na.rm = T), 
                                  se = sd(Fisi, na.rm = T) / sqrt(N.iter),
                                  lHPDI95 = HPDinterval(as.mcmc(Fisi), prob = 0.95, na.rm = T)[1],
                                  uHPDI95 = HPDinterval(as.mcmc(Fisi), prob = 0.95, na.rm = T)[2]))
          
            row.names(outGen) <- rep(NULL, nrow(outGen))
          }
        
        else {
          HPDI95 = 'Inestimable'
          outGen <- rbind(outGen, 
                          cbind(stat = "Na", 
                                mean = mean(Nai, na.rm = T), 
                                se = sd(Nai, na.rm = T) / sqrt(N.iter),
                                lHPDI95 = HPDI95,
                                uHPDI95 = HPDI95))
          outGen <- rbind(outGen, 
                          cbind(stat = "Ne", 
                                mean = mean(Nei, na.rm = T), 
                                se = sd(Nei, na.rm = T) / sqrt(N.iter),
                                lHPDI95 = HPDI95,
                                uHPDI95 = HPDI95))
          outGen <- rbind(outGen, 
                          cbind(stat = "PropPoly", 
                                mean = mean(PropPolyi, na.rm = T), 
                                se = sd(PropPolyi, na.rm = T) / sqrt(N.iter),
                                lHPDI95 = HPDI95,
                                uHPDI95 = HPDI95))
          outGen <- rbind(outGen, 
                          cbind(stat = "He", 
                                mean = mean(Hei, na.rm = T), 
                                se = sd(Hei, na.rm = T) / sqrt(N.iter),
                                lHPDI95 = HPDI95,
                                uHPDI95 = HPDI95))
          outGen <- rbind(outGen, 
                          cbind(stat = "Ho", 
                                mean = mean(Hoi, na.rm = T), 
                                se = sd(Hoi, na.rm = T) / sqrt(N.iter),
                                lHPDI95 = HPDI95,
                                uHPDI95 = HPDI95))
          outGen <- rbind(outGen, 
                          cbind(stat = "IR", 
                                mean = mean(IRi, na.rm = T), 
                                se = sd(IRi, na.rm = T) / sqrt(N.iter),
                                lHPDI95 = HPDI95,
                                uHPDI95 = HPDI95))
          outGen <- rbind(outGen, 
                          cbind(stat = "Fis", 
                                mean = mean(Fisi, na.rm = T), 
                                se = sd(Fisi, na.rm = T) / sqrt(N.iter),
                                lHPDI95 = HPDI95,
                                uHPDI95 = HPDI95))
          
          row.names(outGen) <- rep(NULL, nrow(outGen))      
        }
        
        out <- list(DateTime = self$Date, CompTime = self$SimTime, N.iter = N.iter, N.years = N.years,
                    Lambda = data.frame(mean = c(EmpLambda_mean, exp(StochLogLambda_mean)), 
                                        median = c(EmpLambda_median, exp(StochLogLambda_median)),
                                        lHPDI95 = c(EmpLambda_quant[1], exp(StochLogLambda_quant[1])), 
                                        uHPDI95 = c(EmpLambda_quant[2], exp(StochLogLambda_quant[2])),
                                        l95_se = c(EmpLambda_mean - 1.96*EmpLambda_se, exp(StochLogLambda_mean - 1.96*StochLogLambda_se)), 
                                        u95_se = c(EmpLambda_mean + 1.96*EmpLambda_se, exp(StochLogLambda_mean + 1.96*StochLogLambda_se)),
                                        row.names = c("Empirical Lambda", "Stochastic Lambda")),
                    ExtinctionProb = Prob.extinct,
                    ExtinctionTime = eTime,
                    Pop.size = outSize,
                    Immigrants = imm,
                    GeneticComposition = outGen)
        }
      }
      else {
        out <- list(DateTime = self$Date, CompTime = self$SimTime, N.iter = N.iter, N.years = N.years,
                    Lambda = data.frame(mean = c(EmpLambda_mean, exp(StochLogLambda_mean)), 
                                        median = c(EmpLambda_median, exp(StochLogLambda_median)),
                                        lHPDI95 = c(EmpLambda_quant[1], exp(StochLogLambda_quant[1])), 
                                        uHPDI95 = c(EmpLambda_quant[2], exp(StochLogLambda_quant[2])),
                                        l95_se = c(EmpLambda_mean - 1.96*EmpLambda_se, exp(StochLogLambda_mean - 1.96*StochLogLambda_se)), 
                                        u95_se = c(EmpLambda_mean + 1.96*EmpLambda_se, exp(StochLogLambda_mean + 1.96*StochLogLambda_se)),
                                        row.names = c("Empirical Lambda", "Stochastic Lambda")),
                    ExtinctionProb = Prob.extinct,
                    Pop.size = outSize,
                    Immigrants = imm,
                    GeneticComposition = "No genetic output generated: Check if genOutput = TRUE or if ExtinctionProb = 1")    
      }
      
      out
    },
    
    pullGenoSummary = function(years, genoMetric) {
      Y <- paste("Y", years, sep="")
      o <- list()
      for (y in Y) {
        oGen <- data.frame()
        for (g in genoMetric) {
          if (g == 'PropPoly' | g == 'Ne') {
            stat <- self[[g]][, y]
            oGen <- rbind(oGen,
                          cbind(GenoMetric = g, 
                                mean = mean(stat, na.rm = T), 
                                se = sd(stat, na.rm = T) / sqrt(length(stat)),
                                lHPDI95 = HPDinterval(as.mcmc(stat), prob = 0.95, na.rm = T)[1],
                                uHPDI95 = HPDinterval(as.mcmc(stat), prob = 0.95, na.rm = T)[2]))
          }
          
          else {
            stat <- self[[g]]$mean[, y]
            oGen <- rbind(oGen,
                          cbind(GenoMetric = g, 
                                mean = mean(stat, na.rm = T), 
                                se = sd(stat, na.rm = T) / sqrt(length(stat)),
                                lHPDI95 = HPDinterval(as.mcmc(stat), prob = 0.95, na.rm = T)[1],
                                uHPDI95 = HPDinterval(as.mcmc(stat), prob = 0.95, na.rm = T)[2]))
          }
        }
        o[[y]] <- oGen
      }
      return(o)
    },
    
    immigrants = function () {
      y <- self$years
      ps <- self$populations
      
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
    },
    
    plot = function(fieldStat) {
      if (is.null(fieldStat)) fieldStat <- c('pop.size', 'PropPoly', 'Na', 'Ne', 'He', 'Ho', 'IR', 'Fis', 'lambda', 'extinctTime')
      
      for (p in 1:length(fieldStat)) {
        par(ask=TRUE)
        
        if (fieldStat[p]=='pop.size') {
          ps <- self$pop.size
          for (t in 1:length(ps)) {
            mplots <- list()
            for (i in 1:length(ps[[t]])) {
              psi_mean <- apply(ps[[t]][[i]], 2, function(x) mean(x, na.rm=T))
              psi_low <- apply(ps[[t]][[i]], 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95, na.rm=T)[1])
              psi_hi <- apply(ps[[t]][[i]], 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95, na.rm=T)[2])
              dat_psi <-  data.frame(month = 0:(ncol(ps[[t]][[i]])-1), psi_mean = psi_mean, psi_l95 = psi_low, psi_u95 = psi_hi)
              erib <- aes(ymax = psi_u95, ymin = psi_l95)
              assign(paste('ps_', names(ps[[t]])[i], sep=""), ggplot(dat_psi, aes(x=month, y=psi_mean)) + geom_line(size=1.05) + geom_ribbon(erib, alpha=0.5) +
                       ggtitle(paste(names(ps)[t], 'Population', sep=' ')) +
                       labs(x="Month", y=paste("Population Size:", names(ps[[t]])[i])) + #ylim(c(0,20)) + 
                       theme(axis.text.x=element_text(angle=50, size=10, vjust=0.5),
                             axis.text.y=element_text(size=10),
                             axis.title.x = element_text(size=10, vjust=-0.65),
                             axis.title.y = element_text(size=10, vjust=1)) 
              )
              mplots[[i]] <- get(paste('ps_', names(ps[[t]])[i], sep=""))
            }
            multiplot(plotlist=mplots, cols = 2)
          }
        }
        
        if (fieldStat[p]=='lambda') {
          # By Emp
          lp_Empmean <- apply(self$lambda[, -1], 1, function(x) gm_mean(x, na.rm=T))
          dat_lp <- data.frame(x = lp_Empmean)
          lp1 <- ggplot(data = dat_lp, mapping = aes(x = x)) + geom_histogram(binwidth = 0.0025) + #xlim(c(0.7, 1.1)) +
            aes(y = ..density..) + labs(x="Emp Lambda", y="Density") +
            theme(axis.title.x = element_text(size = 20, vjust = -0.65),
                  axis.title.y = element_text(size = 20, vjust = 1))
          
          # By Last Year
          lp_Stochmean <- exp(apply(self$lambda[, -1], 1, function(x) mean(log(x), na.rm=T)))
          dat_lyp <- data.frame(x = lp_Stochmean)
          lp2 <- ggplot(data = dat_lyp, mapping = aes(x = x)) + geom_histogram(binwidth = 0.0025) + #xlim(c(0.7, 1.1)) +
            aes(y = ..density..) + labs(x="Stoch Lambda", y="Density") +
            theme(axis.title.x = element_text(size = 20, vjust = -0.65),
                  axis.title.y = element_text(size = 20, vjust = 1))
          multiplot(lp1, lp2, cols = 1)
          #par(ask=T) 
        }
        
        if (fieldStat[p]=='extinctTime' & length(na.omit(self$extinctTime)) >= 1) {
          dat_et <- data.frame(x = self$extinctTime)
          etp <- ggplot(data = dat_et, mapping = aes(x = x)) + geom_histogram() + #xlim(c(0.7, 1.1)) +
            aes(y = ..density..) + labs(x="Time to Extinction (Years)", y="Density") +
            theme(axis.title.x = element_text(size = 20, vjust = -0.65),
                  axis.title.y = element_text(size = 20, vjust = 1))
          multiplot(etp)
        }
        
        if (fieldStat[p]=='PropPoly') {
          pp <- self$PropPoly
          pp_mean <- apply(pp, 2, function(x) mean(x, na.rm=T))
          pp_low <- apply(pp, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95, na.rm=T)[1])
          pp_hi <- apply(pp, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95, na.rm=T)[2])
          dat_pp <-  data.frame(year = 0:(length(pp_mean)-1), pp_mean = pp_mean, pp_l95 = pp_low, pp_u95 = pp_hi)
          erib <- aes(ymax = pp_u95, ymin = pp_l95)
          #pp_se <- apply(pp, 2, function(x) sd(x, na.rm=T) / sqrt(nrow(pp)))
          #dat_pp <- data.frame(year = 0:(length(pp_mean)-1), pp_mean = pp_mean, pp_se = pp_se)
          #erib <- aes(ymax = pp_mean + pp_se, ymin = pp_mean - pp_se)
          pp1 <- ggplot(dat_pp, aes(x=year, y=pp_mean)) + geom_line(size=1.05) + geom_ribbon(erib, alpha=0.5) +
            labs(x="Year", y='Prop of Polymorphic Loci') + ylim(c(0,1)) + 
            theme(axis.text.x=element_text(angle=50, size=20, vjust=0.5),
                  axis.text.y=element_text(size=20),
                  axis.title.x = element_text(size=20, vjust=-0.65),
                  axis.title.y = element_text(size=20, vjust=1)) 
          multiplot(pp1, cols=1)
        }
        
        if (fieldStat[p]=='Ne') {
          pp <- self$Ne
          pp_mean <- apply(pp, 2, function(x) mean(x, na.rm=T))
          pp_low <- apply(pp, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95, na.rm=T)[1])
          pp_hi <- apply(pp, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95, na.rm=T)[2])
          dat_pp <-  data.frame(year = 0:(length(pp_mean)-1), pp_mean = pp_mean, pp_l95 = pp_low, pp_u95 = pp_hi)
          erib <- aes(ymax = pp_u95, ymin = pp_l95)
          #pp_se <- apply(pp, 2, function(x) sd(x, na.rm=T) / sqrt(nrow(pp)))
          #dat_pp <- data.frame(year = 0:(length(pp_mean)-1), pp_mean = pp_mean, pp_se = pp_se)
          #erib <- aes(ymax = pp_mean + pp_se, ymin = pp_mean - pp_se)
          pp1 <- ggplot(dat_pp, aes(x=year, y=pp_mean)) + geom_line(size=1.05) + geom_ribbon(erib, alpha=0.5) +
            labs(x="Year", y='Ne') + #ylim(c(0,1)) + 
            theme(axis.text.x=element_text(angle=50, size=20, vjust=0.5),
                  axis.text.y=element_text(size=20),
                  axis.title.x = element_text(size=20, vjust=-0.65),
                  axis.title.y = element_text(size=20, vjust=1)) 
          multiplot(pp1, cols=1)
        }
        
        if (fieldStat[p]=='Na' | fieldStat[p]=='He' | fieldStat[p]=='Ho' | fieldStat[p] =='IR' | fieldStat[p]=='Fis') {
          fi <- self[[fieldStat[p]]]$mean
          fi_mean <- apply(fi, 2, function(x) mean(x, na.rm=T))
          fi_low <- apply(fi, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95, na.rm=T)[1])
          fi_hi <- apply(fi, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95, na.rm=T)[2])
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
    }
))

popClass <- R6Class('popClass',
  public = list(
    popID = NA,
    indsAll = NA,
    indsAlive = NA,
    activeLitters = NA,
    activePairs = NA,
    time = 0,
    pop.size = NA,
    lambda = NA,
    extinct = FALSE,
    Na = NA,
    Ne = NA,
    PropPoly = NA,
    He = NA,
    Ho = NA,
    IR = NA,
    Fis = NA,
    
    initialize = function(popID, ...) {
      self$popID <- popID
    },
    
    startPop = function(startValues, ID, sex, age, mother, father, ageTrans, socialStat, reproStat, genoCols, genOutput) {
      sv <- startValues
      #field('extinct', FALSE)
      self$pop.size <- list(Females = data.frame(M0 = c(0,0,0,0), row.names=c("Kittens", "SubAdults", "Adults", "Total")),
                             Males = data.frame(M0 = c(0,0,0,0), row.names=c("Kittens", "SubAdults", "Adults", "Total")),
                             All = data.frame(M0 = c(0,0,0,0), row.names=c("Kittens", "SubAdults", "Adults", "Total")))
      self$Na <- data.frame(Y0 = c(0,0), row.names=c("mean", "SE"))
      self$Ne <- data.frame(Y0 = c(0))
      self$He <- data.frame(Y0 = c(0,0), row.names=c("mean", "SE"))
      self$Ho <- data.frame(Y0 = c(0,0), row.names=c("mean", "SE"))
      self$IR <- data.frame(Y0 = c(0,0), row.names=c("mean", "SE"))
      self$Fis <- data.frame(Y0 = c(0,0), row.names=c("mean", "SE"))
      self$lambda <- data.frame(Y0 = c(0))
      self$PropPoly <- data.frame(Y0 = c(0))
      self$indsAll <- list()
      self$indsAlive <- list()
      
      # Instantiate individuals
      for (r in 1:nrow(sv)) {
        ind <- indClass$new(animID=sv[r,ID], sex=sv[r,sex], age=sv[r,age], mother=sv[r,mother], father=sv[r,father], socialStat=sv[r,socialStat], 
                            reproStat=sv[r,reproStat], reproHist=as.character(NA), liveStat=TRUE, birthMon=as.numeric(NA), mortMon=as.numeric(NA), 
                            genotype=sv[r,genoCols], censored = F, immigrant = F, ageToAdult=as.numeric(NA))
        self$addToPop(ind)
      }
      self$pullAlive()
      
      # Identify active litters
      iAlive <- self$indsAlive
      aLit <- list() 
      
      # pull subAdults
      sal <- llply(iAlive, function (x) if (x$socialStat == 'SubAdult') x)
      sal <- sal[!sapply(sal, is.null)]
      ageTA <- ageTrans[ageTrans$socialStat == "SubAdult", c("sex", "low", "high")]
      for (sa in sal) {
        if (sa$sex == 'F') {
          rang <- ageTA[ageTA$sex=='F', 'low']:ageTA[ageTA$sex=='F', 'high']
          if (length(rang) == 1) sa$ageToAdult <- rang
          else sa$ageToAdult <- sample(rang, size = 1)
        }
        
        if (sa$sex == 'M') {
          rang <- ageTA[ageTA$sex=='M', 'low']:ageTA[ageTA$sex=='M', 'high']
          if (length(rang) == 1) sa$ageToAdult <- rang
          else sa$ageToAdult <- sample(rang, size = 1)
        }
      }
      
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
      
      self$activeLitters <- aLit
      
      # Start breeding males list
      aMales <- llply(iAlive, function (x) if (x$socialStat == 'Adult' & x$sex == 'M') x)
      aMales <- aMales[!sapply(aMales, is.null)]
      nPairs <- list()
      for (aM in aMales) {
        Male <- aM
        Females <- list()
        newPair <- list(Male = Male, Females = Females)
        nPairs <- append(nPairs, list(newPair))
      }
      self$activePairs <- nPairs
      
      # Update stats
      self$updateStats(genOutput)
    },
    
    addToPop = function(individual) {
      self$indsAll <- append(self$indsAll, individual)
    },
    
    pullAlive = function() {
      alive <- llply(self$indsAll, function(x) if (x$liveStat==TRUE) x)
      alive <- alive[!sapply(alive, is.null)]  
      self$indsAlive <- alive
    },
    
    tabAll = function() {
      dat <- self$indsAll
      out <- c()
      for (r in 1:length(dat)) {
        out = rbind(out, dat[[r]]$tab())
      }
      out
    },
    
    tabAlive = function() {
      dat <- self$indsAlive
      if (length(dat)==0) stop('No individuals are listed as alive.')
      out <- c()
      for (r in 1:length(dat)) {
        out = rbind(out, dat[[r]]$tab())
      }
      out
    },
    
    updateStats = function (genOutput) {
      # pull living individuals
      iAlive <- self$indsAlive
    
      # update population size
      op <- list()
      of <- matrix(0, nrow = 4, ncol=1, byrow = F)
      of[1,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='Kitten' & x$sex=='F'))))
      of[2,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='SubAdult' & x$sex=='F'))))
      of[3,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='Adult' & x$sex=='F'))))
      of[4,1] <- sum(of[1:3,1])
      op$Females <- of
    
      om <- matrix(0, nrow = 4, ncol=1, byrow = F)
      om[1,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='Kitten' & x$sex=='M'))))
      om[2,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='SubAdult' & x$sex=='M'))))
      om[3,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='Adult' & x$sex=='M'))))
      om[4,1] <- sum(om[1:3,1])
      op$Males <- om
    
      ot <- matrix(0, nrow = 4, ncol=1, byrow = F)
      ot[1,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='Kitten'))))
      ot[2,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='SubAdult'))))
      ot[3,1] <- sum(unlist(llply(iAlive, function(x) sum(x$socialStat=='Adult'))))
      ot[4,1] <- sum(ot[1:3,1])
      op$All <- ot
    
      for (subP in 1:length(self$pop.size)) {
        self$pop.size[[subP]][, paste("M", self$time, sep="")] <- op[[subP]]
      }
    
      # update extinction
      if (op$All[4, 1] <= 1) self$extinct <- TRUE
      else {
        if (sum(unlist(llply(iAlive, function(x) sum(x$sex=='M')))) < 1) self$extinct <- TRUE
        if (sum(unlist(llply(iAlive, function(x) sum(x$sex=='F')))) < 1) self$extinct <- TRUE
      }
    
      if ((self$time/12)%%1 == 0) {
        # Assign year for field names
        year <- paste("Y", self$time/12, sep = "")
        
        # lambda
        if (self$time == 0) self$lambda[, year] <- 1
        else {self$lambda[, year] <- self$pop.size$All[nrow(self$pop.size$All), ncol(self$pop.size$All)] / self$pop.size$All[nrow(self$pop.size$All), ncol(self$pop.size$All) - 12]}
      
        if (genOutput) {
          ### Genetic metrics
          g <- pullGenos(iAlive)
          g_genind <- df2genind(createGenInput(g), sep="_", ind.names = g[,1])
        
          if (Sys.info()[[1]] == 'Windows') sink('NUL')
          else {sink('aux')}
          sumGind <- summary(g_genind)
          g_genpop <- genind2genpop(g_genind, pop = rep(self$popID, nrow(g)))
          sink(NULL)
        
          # Allelic richness Na
          self$Na[, year] <- c(mean(g_genind@loc.n.all), sd(g_genind@loc.n.all) / sqrt(length(g_genind@loc.n.all)))
      
          # Ne
          if (year == 'Y0') {
            fem <- llply(iAlive, function (x) if (x$sex == 'F' & x$socialStat == 'Adult') x)
            fem <- length(fem[!sapply(fem, is.null)])
            mal <- llply(iAlive, function (x) if (x$sex == 'M' & x$socialStat == 'Adult') x)
            mal <- length(mal[!sapply(mal, is.null)]) 
            self$Ne[, year] <- (4 * fem * mal) / (fem + mal)
          }
          else {
            fem <- llply(iAlive, function (x) if (x$sex == 'F' & grepl(':', x$reproHist) >= 1) x)
            fem <- length(fem[!sapply(fem, is.null)])
            mal <- llply(iAlive, function (x) if (x$sex == 'M' & grepl(':', x$reproHist) >= 1) x)
            mal <- length(mal[!sapply(mal, is.null)]) 
            self$Ne[, year] <- (4 * fem * mal) / (fem + mal)
          }
        
          # Proportion of alleles polymorphic
          self$PropPoly[, year] <- mean(isPoly(g_genind, by=c("locus")))
        
          # Expected heterozygosity He and Ho
          Hexp <- sumGind$Hexp
          Hobs <- sumGind$Hobs
          self$He[, year] <- c(mean(Hexp), sd(Hexp) / sqrt(length(Hexp)))
          self$Ho[, year] <- c(mean(Hobs), sd(Hobs) / sqrt(length(Hobs)))
        
          # Individual heterozygosity IR
          self$IR[, year] <- c(mean(ir(g[,-1])), sd(ir(g[,-1])) / sqrt(nrow(g)))
        
          # Inbreeding coefficient Fis
          Fis_ind <- sapply(inbreeding(g_genind, N=50), mean)
          self$Fis[, year] <- c(mean(Fis_ind), sd(Fis_ind) / sqrt(length(Fis_ind)))
        }
      }
    },
    
    addImmigrants = function(iP, immMaleProb, ageTrans) {
      ro <- sample(1:nrow(iP), size = 1)
      
      # Pull geno and create new ID
      newgeno <- iP[ro, -c(1:2)]
      newID <- paste('iid', (length(self$indsAll) + 1), sep="")
      
      # Determine sex and age of immigrant
      if (runif(1) <= immMaleProb) {
        sex <- 'M'
        age <- round(runif(1, min=ageTrans[ageTrans$sex == 'M' & ageTrans$socialStat == 'Kitten', 'age'],
                           max=ageTrans[ageTrans$sex == 'M' & ageTrans$socialStat == 'SubAdult', 'age']))
        rang <- ageTrans[ageTrans$sex=='M' & ageTrans$socialStat == "SubAdult", 'low']:ageTrans[ageTrans$sex=='M' & ageTrans$socialStat == "SubAdult", 'high']
        if (length(rang) == 1) ata <- rang
        else ata <- sample(rang, size = 1)   
      }
      else {
        sex <- 'F'
        age <- round(runif(1, min=ageTrans[ageTrans$sex == 'F' & ageTrans$socialStat == 'Kitten', 'age'],
                           max=ageTrans[ageTrans$sex == 'F' & ageTrans$socialStat == 'SubAdult', 'age']))
        rang <- ageTrans[ageTrans$sex=='F' & ageTrans$socialStat == "SubAdult", 'low']:ageTrans[ageTrans$sex=='F' & ageTrans$socialStat == "SubAdult", 'high']
        if (length(rang) == 1) ata <- rang
        else ata <- sample(rang, size = 1)
      }
      
      # Create new immigrant and add to population
      imm <- indClass$new(animID=newID, sex=sex, age=age, mother=as.character(NA), father=as.character(NA), 
                          socialStat="SubAdult", reproStat=FALSE, reproHist=as.character(NA), liveStat=TRUE, censored=FALSE,
                          birthMon=self$time - age, mortMon=as.numeric(NA), genotype=newgeno, immigrant=TRUE, ageToAdult=ata)
      self$addToPop(imm)
      
      # Update living populations
      self$pullAlive()
      
      # Return subset of immigrant populations
      return(iP[-ro, ])
    },
    
    stageAdjust = function(ageTrans, Km, Kf, minMaleReproAge) {
      iAlive <- self$indsAlive
      
      kitsAlive <- llply(iAlive, function(x) if (x$socialStat=='Kitten') x)
      kitsAlive <- kitsAlive[!sapply(kitsAlive, is.null)]
      
      tsubAdultFAlive <- llply(iAlive, function(x) if (x$socialStat=='SubAdult' & x$sex == 'F' & 
                                                       x$age > x$ageToAdult) x)
      tsubAdultFAlive <- tsubAdultFAlive[!sapply(tsubAdultFAlive, is.null)]
      
      tsubAdultMAlive <- llply(iAlive, function(x) if (x$socialStat=='SubAdult' & x$sex == 'M' &
                                                       x$age > x$ageToAdult) x)
      tsubAdultMAlive <- tsubAdultMAlive[!sapply(tsubAdultMAlive, is.null)]
      
      adultFemalesAlive <- llply(iAlive, function(x) if (x$socialStat=='Adult' & x$sex == 'F') x)
      adultMalesAlive <- llply(iAlive, function(x) if (x$socialStat=='Adult' & x$sex == 'M') x)
      adultFemalesAlive <- adultFemalesAlive[!sapply(adultFemalesAlive, is.null)]
      adultMalesAlive <- adultMalesAlive[!sapply(adultMalesAlive, is.null)]
      
      if (length(kitsAlive) > 0) {
        for (k in 1:length(kitsAlive)) {
          kind <- kitsAlive[[k]]
          if (kind$sex == 'F') {
            if (kind$age > ageTrans[ageTrans$sex == 'F' & ageTrans$socialStat == 'Kitten', 'age']) {
              kind$socialStat = "SubAdult" 
              rang <- ageTrans[ageTrans$sex=='F' & ageTrans$socialStat == "SubAdult", 'low']:ageTrans[ageTrans$sex=='F' & ageTrans$socialStat == "SubAdult", 'high']
              if (length(rang) == 1) kind$ageToAdult <- rang
              else kind$ageToAdult <- sample(rang, size = 1)
            }
          }
          else {
            if (kind$age > ageTrans[ageTrans$sex == 'M' & ageTrans$socialStat == 'Kitten', 'age']) {
              kind$socialStat = "SubAdult" 
              rang <- ageTrans[ageTrans$sex=='M' & ageTrans$socialStat == "SubAdult", 'low']:ageTrans[ageTrans$sex=='M' & ageTrans$socialStat == "SubAdult", 'high']
              if (length(rang) == 1) kind$ageToAdult <- rang
              else kind$ageToAdult <- sample(rang, size = 1)          
            }
          }
        }
      }
      
      if (length(tsubAdultFAlive) > 0) {
        adF <- length(adultFemalesAlive)
        
        # Determine the K for the month
        if (adF < Kf[1, 1]) {
          rown <- which(runif(1) <= cumsum(Kf[1, -1]))
          kf <- Kf[rown[1], 1]
        }    
        else {
          rown <- which(runif(1) <= cumsum(Kf[max(which(Kf[, 1] <= adF)), -1]))
          kf <- Kf[rown[1], 1]
        }
        
        # Identify the allowed space
        allowF <- kf - adF
        
        # If reduction in K, remove adults
        if (allowF < 0) {
          numRemove <- abs(allowF) 
          for (nR in 1:numRemove) {
            ages <- unlist(llply(adultFemalesAlive, function(x) x$age))
            toRemove <- which(ages == min(ages))
            if (length(toRemove) > 1) toRemove <- sample(toRemove, size=1)
            adultFemalesAlive[[toRemove]]$liveStat <- FALSE
            adultFemalesAlive[[toRemove]]$mortMon <- self$time
            adultFemalesAlive <- adultFemalesAlive[[-toRemove]]
          }
        }
        
        # If no space available, kill off subAdult females of the appropriate age
        if (allowF <= 0) {
          invisible(llply(tsubAdultFAlive, function (x) {
            x$liveStat <- FALSE
            x$mortMon <- self$time
            x$censored <- TRUE
          }))
        }
        
        # Else if space is available, sample subAdult females of the appropriate age
        else {
          sampF.size <- min(length(tsubAdultFAlive), allowF)
          samp <- sample(1:length(tsubAdultFAlive), size = sampF.size)
          invisible(llply(tsubAdultFAlive[samp], function(x) {
            x$socialStat = 'Adult'
            x$reproStat = TRUE
          }))
          invisible(llply(tsubAdultFAlive[-samp], function(x) {
            x$liveStat = FALSE
            x$mortMon <- self$time
            x$censored <- TRUE
          }))
        }
      }
      
      if (length(tsubAdultMAlive) > 0) {
        adM <- length(adultMalesAlive)
        
        # Determine the K for the month
        if (adM < Km[1, 1]) {
          rown <- which(runif(1) <= cumsum(Km[1, -1]))
          km <- Km[rown[1], 1]
        }
        else {
          rown <- which(runif(1) <= cumsum(Km[max(which(Km[, 1] <= adM)), -1]))
          km <- Km[rown[1], 1]
        }
        
        # Identify the allowed space
        allowM <- km - adM
        
        # If reduction in K, remove adults
        if (allowM < 0) {
          numRemove <- abs(allowM) 
          for (nR in 1:numRemove) {
            ages <- unlist(llply(adultMalesAlive, function(x) x$age))
            toRemove <- which(ages == min(ages))
            if (length(toRemove) > 1) toRemove <- sample(toRemove, size=1)
            adultMalesAlive[[toRemove]]$liveStat <- FALSE
            adultMalesAlive[[toRemove]]$mortMon <- self$time
            adultMalesAlive <- adultMalesAlive[[-toRemove]]
          }
        }
        
        # If no space available, kill off subAdult males of the appropriate age
        if (allowM <= 0) {
          invisible(llply(tsubAdultMAlive, function (x) {
            x$liveStat <- FALSE
            x$mortMon <- self$time
            x$censored <- TRUE
          }))
        }
        
        # Else if there is space, sample subAdult males of the appropriate age
        else {
          sampM.size <- min(length(tsubAdultMAlive), allowM)
          samp <- sample(1:length(tsubAdultMAlive), size = sampM.size)
          invisible(llply(tsubAdultMAlive[samp], function(x) {
            x$socialStat = 'Adult'
            #x$reproStat = TRUE
          }))
          invisible(llply(tsubAdultMAlive[-samp], function(x) {
            x$liveStat = FALSE
            x$mortMon <- self$time
            x$censored <- TRUE
          }))
        }
      }
      
      #  If no adult or subadult males of the appropriate age, allow younger males to move up
      if (length(tsubAdultMAlive) == 0 & length(adultMalesAlive) == 0) {
        tsaMAlive <- llply(iAlive, function(x) if (x$socialStat=='SubAdult' & x$sex == 'M' &
                                                   x$age >= minMaleReproAge) x)
        tsaMAlive <- tsaMAlive[!sapply(tsaMAlive, is.null)]
        
        if (length(tsaMAlive) > 0) {
          rown <- which(runif(1) <= cumsum(Km[1, -1]))
          allowM <- Km[rown[1], 1]
          sampM.size <- min(length(tsaMAlive), allowM)
          samp <- sample(1:length(tsaMAlive), size = sampM.size)
          invisible(llply(tsaMAlive[samp], function(x) {
            x$socialStat = 'Adult'
            #x$reproStat = TRUE
          }))
        }
      }
      
      self$pullAlive()
    },
    
    updateBreedStat = function(ageTrans, maxN_ReproMale) {
      
      # Update male breed stat
      m_repro <- llply(self$indsAlive, function(x) if (x$sex=="M" & x$socialStat=="Adult" & x$reproStat==TRUE) x)
      m_repro <- m_repro[!sapply(m_repro, is.null)]  
      
      if (length(m_repro) < maxN_ReproMale) {
        m_alive <- llply(self$indsAlive, function(x) if (x$sex=="M" & x$socialStat=="Adult" & x$reproStat==FALSE) x)
        m_alive <- m_alive[!sapply(m_alive, is.null)] 
        
        if (length(m_alive) != 0) {
          nMales <- maxN_ReproMale - length(m_repro)

          if (length(m_alive) > nMales) {
            invisible(llply(sample(m_alive, size = nMales), function(x) {
              x$reproStat = TRUE
            }))
          }
          else {
            invisible(llply(m_alive, function(x) {
              x$reproStat = TRUE
            })) 
          }
        }
      }
      
      ## Females with litters
      aL <- self$activeLitters
      
      # Adjust litters based on mortality of mother
      if (length(aL) > 0 ) {
        for (l in length(aL):1) {
          if (aL[[l]]$mother$liveStat == FALSE) {
            # kill any remaining kittens
            if (length(aL[[l]]$kittens) > 0) {
              for (ks in 1:length(aL[[l]]$kittens)) {
                if (aL[[l]]$kittens[[ks]]$age < 12) {
                  aL[[l]]$kittens[[ks]]$liveStat <- FALSE
                  aL[[l]]$kittens[[ks]]$mortMon <- self$time
                }
                else {
                  aL[[l]]$kittens[[ks]]$socialStat <- 'SubAdult'
                  rang <- ageTrans[ageTrans$sex==aL[[l]]$kittens[[ks]]$sex & ageTrans$socialStat == "SubAdult", 'low']:ageTrans[ageTrans$sex==aL[[l]]$kittens[[ks]]$sex & ageTrans$socialStat == "SubAdult", 'high']
                  if (length(rang) == 1) aL[[l]]$kittens[[ks]]$ageToAdult <- rang
                  else aL[[l]]$kittens[[ks]]$ageToAdult <- sample(rang, size = 1)  
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
        self$activeLitters <- aL
      }
    },
    
    reproduce = function(litterProbs,probBreed,probFemaleKitt,lociNames) {
      # Generate monthly probability of breeding
      tPB <- betaval(probBreed$prob, probBreed$se)
      
      # Pull reproductive adults
      m_alive <- llply(self$indsAlive, function(x) if (x$sex=="M" & x$socialStat=="Adult" & x$reproStat==TRUE) x)
      m_alive <- m_alive[!sapply(m_alive, is.null)]  
      
      ######################
      ### Managing pairs ###
      ######################
      aPairs <- self$activePairs
      
        # Males
      alive <- unlist(llply(1:length(aPairs), function(x) aPairs[[x]]$Male$liveStat))
      new_males_ind <- which(!(m_alive %in% llply(aPairs, function(x) x$Male)))
        
        # Replace dead males
      if (any(!alive)) {
        dead_index <- which(!alive)
        for (di in dead_index) {
          if (any(new_males_ind)) {
            aPairs[[di]]$Male <- m_alive[[new_males_ind[1]]]
            new_males_ind <- new_males_ind[-1]            
          }
          else {aPairs[[di]]$Male <- list()}
        }
      }
      
        # Add new males
      if (any(new_males_ind)) {
        for (nm in new_males_ind) {
          aPairs <- append(aPairs, list(list(Male = m_alive[[nm]], Females = list())))
        }           
      }     
      
        # Females
      #f_alive <- llply(self$indsAlive, function(x) if (x$sex=="F" & x$socialStat=="Adult" & x$reproStat==TRUE) x)
      f_alive <- llply(self$indsAlive, function(x) if (x$sex=="F" & x$socialStat=="Adult") x)
      f_alive <- f_alive[!sapply(f_alive, is.null)] 
      new_females_ind <- which(!(f_alive %in% llply(aPairs, function(x) x$Females)))
      
        # Remove dead females

      
      ## Produce new litters
      #if (length(f_alive) == 0 | length(m_alive) == 0) field('extinct', TRUE)
      if (length(f_alive) != 0 & length(m_alive) != 0) {
        for (f in length(f_alive):1) {
          if (runif(1) <= tPB) {
            # Select male mate
            mate <- sample(m_alive, size = 1)  # For random mating
            
            # Generate number of kitts
            #numKitts <- littSize(litterProbs)
            #if (signif(sum(litterProbs$prob)) != 1) 
            #  stop("Litter size probabilities must sum to 1")
            indProb <- runif(1)
            high <- min(which(litterProbs$cumProbs > indProb))
            numKitts <- litterProbs[high, 'LitterSize']
            
            # Breed
            #f_alive[[f]]$femBreed(mate[[1]], numKitts, probFemaleKitt, lociNames, pop1)
            f_alive[[f]]$femBreed(mate[[1]], numKitts, probFemaleKitt, lociNames, self)
          }
        }
        self$pullAlive()
      }
    },
    
    kill = function(surv) {
      tSurv <- newSurv(surv)
      alive <- self$indsAlive
      
      # Assess survival
      if (length(alive) > 0) {
        for (i in 1:length(alive)) {
          ind <- alive[[i]]
          si <- tSurv[tSurv$sex == ind$sex & tSurv$socialStat == ind$socialStat, 's']
          ind$liveStat <- runif(1) <= si
          
          #update mortMon for individual
          if (ind$liveStat == FALSE) ind$mortMon <- self$time
        }
      }
      
      self$pullAlive()
    },
    
    incremTime = function(senesc) {
      self$time <- self$time + 1
      
      # age individuals
      alive <- self$indsAlive
      if (length(alive) > 0) {
        for (i in 1:length(alive)) {
          alive[[i]]$age <- alive[[i]]$age + 1
          if (alive[[i]]$age > (senesc * 12)) alive[[i]]$liveStat <- FALSE
        }
      }
      self$pullAlive()
    }
))

indClass <- R6Class('indClass',
  public = list(
    animID = NA,
    sex = NA,
    age = NA,
    mother = NA,
    father = NA,
    socialStat = NA,
    ageToAdult = NA,
    reproStat = FALSE,
    reproHist = NA,
    liveStat = TRUE,
    censored = FALSE,
    birthMon = NA,
    mortMon = NA,
    immigrant = FALSE,
    genotype = NA,
    
    ## Methods
    initialize = function(animID, sex, age, mother, father, socialStat, ageToAdult, reproStat, reproHist, liveStat, censored, 
                          birthMon, mortMon, immigrant, genotype) {
      if (!missing(animID)) self$animID <- animID
      if (!missing(sex)) self$sex <- sex
      if (!missing(age)) self$age <- age
      if (!missing(mother)) self$mother <- mother
      if (!missing(father)) self$father <- father
      if (!missing(socialStat)) self$socialStat <- socialStat
      if (!missing(ageToAdult)) self$ageToAdult <- ageToAdult
      if (!missing(reproStat)) self$reproStat <- reproStat
      if (!missing(reproHist)) self$reproHist <- reproHist
      if (!missing(liveStat)) self$liveStat <- liveStat
      if (!missing(censored)) self$censored <- censored
      if (!missing(birthMon)) self$birthMon <- birthMon
      if (!missing(mortMon)) self$mortMon <- mortMon
      if (!missing(immigrant)) self$immigrant <- immigrant
      if (!missing(genotype)) self$genotype <- genotype},
    
    tab = function(...) {
      return(data.frame(animID = self$animID, sex = self$sex, age = self$age, mother = self$mother, father = self$father,
                        socialStat = self$socialStat, ageToAdult = self$ageToAdult, reproStat = self$reproStat,
                        reproHist = self$reproHist, liveStat = self$liveStat, censored = self$censored, birthMon = self$birthMon,
                        mortMon = self$mortMon, immigrant = self$immigrant, self$genotype))
    },
    
    femBreed = function(male, numKittens, probFemaleKitt, lociNames, population) {
      if(self$sex != "F")
        stop("Input mother is not Female")
      if(self$liveStat != TRUE)
        stop("Input mother is not alive and cannot reproduce!")
      if(self$reproStat != TRUE)
        stop("Input mother is not reproductive")
      
      if(male$sex != "M")
        stop("Input father is not male")
      if(male$liveStat != TRUE)
        stop("Input father is not alive and cannot reproduce!")
      if(male$reproStat != TRUE)
        stop("Input father is not reproductive")
      
      pop <- length(population$indsAll)
      
      # Generate new individuals
      # Determine IDs
      idKitt <- paste('sid', seq(pop + 1, pop + numKittens), sep = "")
      
      # Determine sex
      sexKitt <- runif(numKittens, min = 0, max = 1) <= probFemaleKitt
      sexKitt <- ifelse(sexKitt==TRUE, "F", "M")
      
      # Determine genotype...completely random following Mendelian principles
      gts <- rbind(self$genotype, male$genotype)
      
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
        ind <- indClass$new(animID=idKitt[k], sex=sexKitt[k], age=0, mother=self$animID, father=male$animID, socialStat="Kitten", 
                            reproStat=FALSE, reproHist=as.character(NA), liveStat=TRUE, censored=FALSE,
                            birthMon=bm, mortMon=as.numeric(NA), genotype=genoKitt[k,], immigrant=FALSE, ageToAdult=as.numeric(NA))
        
        # add to population
        population$addToPop(ind) 
        kittens <- append(kittens, ind)
      }
      
      # update individuals alive
      #population$pullAlive()
      
      # update activeLitters
      population$activeLitters <- append(population$activeLitters, 
                                         list(list(mother=self, kittens = kittens, gestation = 0)))
      
      # Update reproHist and reproStat for mother and father
      self$reproStat <- FALSE
      self$reproHist <- paste(self$reproHist, numKittens, sep=":")
      male$reproHist <- paste(male$reproHist, numKittens, sep=":")
    }
))