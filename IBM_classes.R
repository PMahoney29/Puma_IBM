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

IDs = newInds(newIDs, sexes, genotypes, ages, socialStats, reproStats, birthMons)
  

### Define Classes

# Individual
setClass("ind",
         representation(genotype = "character",
                        sex = "character",
                        ID = "numeric",
                        age = "numeric",
                        socialStat = "character",
                        reproStat = "logical",
                        birthMon = "numeric",
                        mortMon = "numeric",
                        liveStat = "logical"))

# Class Functions
newInds <- function(newIDs, sexes, genotypes, ages, socialStats, reproStats, birthMons) {
  
  nrecs <- length(newIDs)
  if (nrecs != length(sexes) || nrecs != length(genotypes) || nrecs != length(ages) || nrecs != length(socialStats) 
      || nrecs != length(reproStats) || nrecs != length(birthMons))
    stop("New individual state variables not of equal length")
  
  if (!is.character(genotypes) || !is.character(sexes) || !is.numeric(newIDs) || !is.numeric(ages) 
      || !is.character(socialStats) || !is.logical(reproStats) || !is.numeric(birthMons))
    stop("New individual state variable classes are not properly defined")
  
  mortMon = as.numeric(rep(NA, nrecs))
  liveStat = as.logical(rep("TRUE", nrecs))
  
  new("ind", ID = as.vector(newIDs),
             sex = as.vector(sexes),
             genotype = as.vector(genotypes),
             age = as.vector(ages),
             socialStat = as.vector(socialStats),
             reproStat = as.vector(reproStats),
             birthMon = as.vector(birthMons),
             mortMon = as.vector(mortMon),
             liveStat = as.vector(liveStat)
      )
}


# Accessor functions
printIDs <- function(individuals) individuals@ID
printGeno <- function(individuals) individuals@genotype
printSex <- function(individuals) individuals@sex
printAge <- function(individuals) individuals@age
printSocial <- function(individuals) individuals@socialStat
printRepro <- function(individuals) individuals@reproStat
printBirthMon <- function(individuals) individuals@birthMon
printLive <- function(individuals) individuals@liveStat

# Method functions
setMethod(show, signature("ind"),
          function(object)
            print(data.frame(ID = printIDs(object),
                             Genotype = printGeno(object),
                             Sex = printSex(object),
                             Age = printAge(object),
                             Social = printSocial(object),
                             Repro = printRepro(object),
                             BirthMon = printBirthMon(object),
                             Alive = printLive(object))))

setMethod("Arith", signature(e1 = "ind",
                             e2 = "ind"),
          function(e1, e2)
          {
            #if (!sameloc(e1, e2))
            #  stop("identical locations required")
            newInds(printIDs(e1),
                    printGeno(e1),
                    callGeneric(printAge(e1),
                                printAge(e2)))
          })
