revDiff <- function(current, target) {
  return(current / target - 1)
}

normalize <- function(xFr) {
  return(xFr / sum(xFr))
}

checkElemNames <- function(xFr) {
  if (is.null(names(xFr))) {
    stop("There is no `names()` attribute on the composition vector")
  }
}

compareElemNames <- function(...) {
  arguments <- list(...)
  if (length(arguments) < 2) {
    stop("Less than two component lists are passed to the function")
  }
  result <- lapply(arguments, function(x) all(x == arguments[[1]]))
  if (any(result == FALSE)) {
    stop("The lists of components in the compositions do not match")
  }
}

checkCompStake <- function(xFr) {
  if (any(xFr < 0)) {
    stop("The concentration of one or more components is negative")
  } else if (xFr %>% sum() %>% round(6) == 1) {
    return(xFr)
  } else if (xFr %>% sum() %>% round(4) == 100) {
    message("The % composition is passed to the function,
            it will be returned as a fraction")
    return(xFr / 100)
  } else stop("Passed non-normalized composition, use normalize()")
}

completenessCompCheck <- function(xFr, elemNames) {
  checkElemNames(xFr)
  elemMask <- names(xFr) %in% elemNames
  completeComp <- rep(0, times = length(elemNames))
  names(completeComp) <- elemNames
  completeComp[elemNames %in% names(xFr)] <- xFr[elemMask]
  return(completeComp)
}

molFrToWtFr <- function(molFr, mmElem) {
  return((molFr * mmElem) %>% normalize())
}

wtFrToMolFr <- function(wtFr, mmElem) {
  return((wtFr / mmElem) %>% normalize())
}

additiveProperty <- function(xFr, propElem) {
  xFr <- xFr %>% matrix(nrow = 1)
  propElem <- propElem %>% matrix(ncol = 1)
  property <- xFr %*% propElem %>% as.double()
  return(property)
}

idealGasVM <- function(p, t, R) {
  return(R * t / p)
}

idealGasG <- function(mmSubstance, mmReference) {
  return(mmSubstance / mmReference)
}

idealGasD <- function(mm, p, t, R) {
  return(mm / idealGasVM(p, t, R))
}

realGasZ <- function(molFr, sElem) {
  return(1 - additiveProperty(molFr, sElem) ^ 2)
}

realGasReZ <- function(z, p) {
  return(1 - p / 101.325 * (1 - z))
}

realGasVM <- function(z, p = 101.325, t = 293.15, R = 8.31451) {
  return(z * idealGasVM(p ,t, R))
}

realGasG <- function(mmSubstance, zSubstance, mmReference, zReference) {
  return(idealGasG(mmSubstance, mmReference) * (zReference / zSubstance))
}

realGasD <- function(mm, z, p = 101.325, t = 293.15, R = 8.31451) {
  return(idealGasD(mm, p, t, R) / z)
}

mrcMMrestSystem <- function(gorMass, wtFrRestVapor, mmRestVapor,
                            wtFrRestLiquid, mmRestLiquid) {
  return((wtFrRestVapor * gorMass + wtFrRestLiquid) /
           (wtFrRestVapor * gorMass / mmRestVapor +
              wtFrRestLiquid / mmRestLiquid))
}

mrcMolSystem <- function(gorMass, wtFrVapor, wtFrLiquid, mmElemSystem) {
  return((wtFrVapor * gorMass + wtFrLiquid) / mmElemSystem)
}

mrcMolSystemByMM <- function(gorMass, mmVapor, mmLiquid) {
  return(1 / mmLiquid + gorMass / mmVapor)
}

# The `checkCompStake()` is executed in case of an incorrect gorMass value
mrcWtFrSystem <- function(gorMass, wtFrVapor, wtFrLiquid) {
  return(((wtFrVapor * gorMass + wtFrLiquid) / (1 + gorMass)) %>%
           checkCompStake())
}


mrcMolFrSystem <- function(molElemSystem, molSystem=sum(molElemSystem)) {
  return(molElemSystem / molSystem)
}

mrcMMbyMolSystem <- function(gorMass, molSystem) {
  return((1 + gorMass) / molSystem)
}

reMRcWtFrVapor <- function(gorMass, wtFrSystem, wtFrLiquid) {
  return(wtFrSystem + (wtFrSystem - wtFrLiquid) / gorMass)
}
