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

########
## Test data
########
datum <- read.csv('./Data//genotypes//TestGenData.csv', skip = 2)
geninput <- datum[,1:3]
cols <- seq(ncol(geninput)+1, 110, by=2)
for (c in cols) {
 geninput <- cbind(geninput,paste(datum[,c], datum[,c+1],sep="_")) 
 names(geninput)[ncol(geninput)] <- names(datum)[c]
}
genID <- df2genind(geninput[-c(1:3)], sep="_", ind.names=geninput$ID, loc.names=names(geninput[-c(1:3)]), missing = NA, type='codom') 
#mySamp <- genID[sample(1:nrow(genID@tab), 10)]
#mySamp
genout <- genind2df(genID, oneColPerAll=TRUE)
write.csv(genout, "./Data/genotypes/startValues.csv")

# Test instances of the above two classes
startValues <- read.csv('./Data/genotypes/startValues.csv', stringsAsFactors=F)
genoCols = 6:ncol(startValues); ID = "ID"; sex = 'sex'; age = 'age'; socialStat = 'socialStat'; reproStat = 'reproStat';
pop1 <- popClass$new(popID = 'Population_1', time=0)
pop1$startPop(startValues=startValues, ID='ID', sex='sex', age='age', socialStat='socialStat', reproStat='reproStat', genoCols=genoCols)


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
    socialStat = 'character',
    reproStat = 'logical',
    reproHist = 'numeric',
    liveStat = 'logical',
    birthMon = 'numeric',
    mortMon = 'ANY',         ## Need to fix to be either numeric or NA
    genotype = 'data.frame'))


## indClass methods
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
  if(class(popName)[1]!="popClass") 
    stop("Population object must be of class : 'popClass'")

  ## need to test for unique names
  
  # extending the list of individuals
  popName$indsAll <- append(popName$indsAll, .self)
})


## popClass methods
popClass$methods(startPop = function(startValues, ID, sex, age, socialStat, reproStat, genoCols) {
  sv <- startValues
  for (r in 1:nrow(sv)) {
    ind <- indClass$new(animID=sv[r,ID], sex=sv[r,sex], age=sv[r,age], socialStat=sv[r,socialStat], reproStat=sv[r,reproStat], reproHist=0, 
                        liveStat=TRUE, birthMon=0, mortMon=NA, genotype=sv[r,genoCols])
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
  })

  # Update lambda






