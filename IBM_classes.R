############################################################################
#  Individually based model investigating the genetic consequences of small
#  populations.  Started - 2/2015
#  Author: Peter Mahoney, USU PhD Student
############################################################################


########
## Test data
########
newIDs = 1:5
sexes = c("F","M","M","F","F")
genotypes = c("CCaaBB","CCaaBB","CCaaBB","CCaaBB","CCaaBB")
ages = rep(0,5)
socialStats = rep("Pup",5)
reproStats = rep(FALSE, 5)
birthMons = rep(0,5)
mortMons = rep(-999,5)
liveStats = rep(TRUE,5)

IDs = newInds(newIDs, sexes, genotypes, ages, socialStats, reproStats, birthMons, mortMons, liveStats)
  

### Define Classes

# Individual
setClass("individuals",
         representation(Genotype = "character",
                        Sex = "character",
                        ID = "numeric",
                        Age = "numeric",
                        SocialStat = "character",
                        ReproStat = "logical",
                        BirthMon = "numeric",
                        MortMon = "numeric",
                        LiveStat = "logical"))

# Class Functions
newInds <- function(newIDs, sexes, genotypes, ages, socialStats, reproStats, birthMons, mortMons, liveStats) {
  
  nrecs <- length(newIDs)
  if (nrecs != length(sexes) || nrecs != length(genotypes) || nrecs != length(ages) || nrecs != length(socialStats) 
      || nrecs != length(reproStats) || nrecs != length(birthMons) || nrecs != length(mortMons) || nrecs != length(liveStats))
    stop("New individual state variables not of equal length")
  
  if (!is.character(genotypes) || !is.character(sexes) || !is.numeric(newIDs) || !is.numeric(ages) 
      || !is.character(socialStats) || !is.logical(reproStats) || !is.numeric(birthMons) || !is.numeric(mortMons) || !is.logical(liveStats))
    stop("New individual state variable classes are not properly defined")
  
  #mortMon = as.numeric(rep(NA, nrecs))
  #liveStat = as.logical(rep("TRUE", nrecs))
  
  new("individuals", ID = as.vector(newIDs),
             Sex = as.vector(sexes),
             Genotype = as.vector(genotypes),
             Age = as.vector(ages),
             SocialStat = as.vector(socialStats),
             ReproStat = as.vector(reproStats),
             BirthMon = as.vector(birthMons),
             MortMon = as.vector(mortMons),
             LiveStat = as.vector(liveStats)
      )
}


# Accessor functions
printIDs <- function(individuals) individuals@ID
printGeno <- function(individuals) individuals@Genotype
printSex <- function(individuals) individuals@Sex
printAge <- function(individuals) individuals@Age
printSocial <- function(individuals) individuals@SocialStat
printRepro <- function(individuals) individuals@ReproStat
printBirthMon <- function(individuals) individuals@BirthMon
printMortMon <- function(individuals) individuals@MortMon
printLive <- function(individuals) individuals@LiveStat

# Method functions
setMethod(show, signature("individuals"),
          function(object)
            print(data.frame(ID = printIDs(object),
                             Genotype = printGeno(object),
                             Sex = printSex(object),
                             Age = printAge(object),
                             Social = printSocial(object),
                             Repro = printRepro(object),
                             BirthMon = printBirthMon(object),
                             MortMon = printMortMon(object),
                             Alive = printLive(object))))

setMethod("[",
          signature(x = "individuals",
                    i = "ANY",
                    j = "missing",
                    drop = "missing"),
          function(x, i, j)
            newInds(printIDs(x)[i],
                    printGeno(x)[i],
                    printSex(x)[i],
                    printAge(x)[i],
                    printSocial(x)[i],
                    printRepro(x)[i],
                    printBirthMon(x)[i],
                    printMortMon(x)[i],
                    printLive(x)[i]))

setGeneric("addInds",
           function(a, b)
             standardGeneric("addInds")
           )

setMethod("addInds",
          signature("individuals", "individuals"),
          function(a, b) {
            names <- slotNames(a)
            for (name in names) {
              slot(a, name) <- c(slot(a, name), slot(b, name))
            }
            if (anyDuplicated(a@ID)) stop("Non-unique animal IDs created")
            return(a)
          })

## IDs[printLive(IDs)==TRUE]

setGeneric("Arith", signature(e1 = "individuals",
                             e2 = "individuals"),
          function(e1, e2)
          {
            #if (!sameloc(e1, e2))
            #  stop("identical locations required")
            newInds(printIDs(e1),
                    printGeno(e1),
                    callGeneric(printAge(e1),
                                printAge(e2)))
          })
