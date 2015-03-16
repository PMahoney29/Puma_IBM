############################################################################
#  Individually based model classes
#  investigating the genetic consequences of small populations.  
#  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
############################################################################

#Load required packages
library(methods)
library(adegenet)
library(plyr)


########
## Test data
########
animID = 'a1'
sex = "F"
age = 0
socialStat = 'Pup'
reproStat = FALSE
reproHist = 0
birthMon = 0
mortMon = as.numeric(NA)
liveStat = TRUE
genotype = data.frame(L1="Aa", L2="cD", stringsAsFactors=FALSE)


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
    months = 'numeric',
    pop.size = 'numeric',
    lambda = 'numeric'
  ))

indClass <- setRefClass(
  Class = 'indClass',
  fields = list(
    animID = 'character',
    sex = 'character',
    age = 'numeric',
    socialStat = 'character',
    reproStat = 'logical',
    reproHist = 'numeric',
    liveStat = 'logical',
    birthMon = 'numeric',
    mortMon = 'numeric',
    genotype = 'data.frame'))

#,
#  prototypes = list(
#    age = 0,
#    socialStat = "Pup",
#    reproStat = FALSE,
#    reproHist = 0,
#    liveStat = TRUE,
#    birthMon = 0,
#    mortMon = NA)
#  )

## indClass methods
  # Test instances of the above two classes
a1 <- indClass$new(animID=animID, sex=sex, age=age, socialStat=socialStat, 
                   reproStat=reproStat, reproHist=reproHist, liveStat=liveStat, birthMon=birthMon, mortMon=mortMon, genotype=genotype)
a2 <- indClass$new(animID="a2", sex=sex, age=age, socialStat=socialStat, 
                   reproStat=reproStat, reproHist=reproHist, liveStat=liveStat, birthMon=birthMon, mortMon=mortMon, genotype=genotype)
a3 <- indClass$new(animID="a3", sex=sex, age=age, socialStat=socialStat, 
                   reproStat=reproStat, reproHist=reproHist, liveStat=liveStat, birthMon=birthMon, mortMon=mortMon, genotype=genotype)
pop1 <- popClass$new(popID = 'Population_1')
a1$addToPop(pop1)
a2$addToPop(pop1)
a3$addToPop(pop1)
pop1$tabIndsAll()
pop1$pullAlive()
pop1$tabAlive()

  # Print method for individual data
indClass$methods(tab = function() {
  fields <- names(.refClassDef@fieldClasses)
  out <- data.frame()
  for (fi in fields) {
    #####  Will need to fix depending on what I decide to do with genind() object classes
    if (class(field(fi))=='data.frame') {
      g <- c(field(fi), sep = " ")
      out[1,fi] <- do.call(paste, g)
    }
    else {out[1,fi] <- field(fi)}
  }
  names(out) <- fields
  out
 })


  # Add method for including individual data to popClass object
indClass$methods(addToPop = function(popName) {
  if(class(popName)[1]!="popClass") 
    stop("Population object must be of class : 'popClass'")

  ## need to test for unique names
  
  # extending the list of individuals
  popName$indsAll <- append(popName$indsAll, .self)
  
  # updating population size
  #popName$pop.size <- length(popName$inds)
})


## popClass methods
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

  # Update months
  
  # Update lambda

## aliveClass methods


