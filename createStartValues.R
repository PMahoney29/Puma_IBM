######################
######################
##  Generating the starting population file for CalPuma IBM PVA
######################
######################

#Load required packages
library(methods)
library(adegenet)
library(plyr)
source('./classes_IBM.R')

########
## Initial Values modification
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

startValues <- read.csv('./Data/genotypes/startValues.csv', stringsAsFactors=F)
lociNames <- unique(sub("[.].*$","",names(startValues)[-c(1:5)]))

genoCols = 6:ncol(startValues); startValues$age <- as.numeric(startValues$age); 
pop1 <- popClass$new(popID = 'Population_1', time=0)
pop1$startPop(startValues=startValues, ID='ID', sex='sex', age='age', socialStat='socialStat', reproStat='reproStat', genoCols=genoCols)

testF23 <- pop1$indsAlive[[6]]
testF13 <- pop1$indsAlive[[3]]
testF19 <- pop1$indsAlive[[4]]
testM12 <- pop1$indsAlive[[2]]
testM22 <- pop1$indsAlive[[5]]

testF23$femBreed(testM12, 1, 1, lociNames, pop1)
testF23$femBreed(testM22, 1, 0, lociNames, pop1)
testF13$femBreed(testM22, 1, 1, lociNames, pop1)
testF13$femBreed(testM12, 1, 0, lociNames, pop1)

testF19$femBreed(testM12, 1, 1, lociNames, pop1)
pop1$indsAll[[14]]$socialStat <- 'SubAdult'
testF13$femBreed(testM22, 1, 0, lociNames, pop1)
pop1$indsAll[[15]]$socialStat <- 'SubAdult'

out <- pop1$tabAlive()
write.csv(out, "startValues_complete.csv", row.names=FALSE, na = "")
##########
##########