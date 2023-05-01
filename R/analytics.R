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
  } else {
    stop("Passed non-normalized composition, use normalize()")
  }
}

completenessCompCheck <- function(xFr, elemNames) {
  checkElemNames(xFr)
  elemMask <- names(xFr) %in% elemNames
  completeComp <- rep(0, times = length(elemNames))
  names(completeComp) <- elemNames
  completeComp[elemNames %in% names(xFr)] <- xFr[elemMask]
  return(completeComp)
}

checkRealGasZ <- function(z) {
  if (!((z > 0.9) && (z <= 1))) {
    warning("z-factor outside the valid value range: 0.9 < z <= 1")
  }
}

checkRealGasP <- function(p) {
  if (!((p >= 90) && (p <= 110))) {
    warning("Pressure outside the permissible value range: 90 <= p <= 110")
  }
}

checkRealGasT <- function(t) {
  tAllowable <- c(273.15, 288.15, 288.70, 293.15)
  if (!(t %in% tAllowable)) {
    warning("Temperature is not equal to standard: 273.15, 288.15, 288.70, 293.15 k")
  }
}

reMRcSumVapor <- function(compSetVapor, gorVol, wtFrSystem,
                          wtFrLiquid, dLiquid) {
  return(function(dVapor) {
    checkElemNames(wtFrSystem)
    wtFrSystem %<>% checkCompStake()
    checkElemNames(wtFrLiquid)
    wtFrLiquid %<>% checkCompStake()
    compareElemNames(names(wtFrSystem), names(wtFrLiquid))

    return(reMRcWtFrVapor(wtFrSystem, wtFrLiquid,
                          gorVol * dVapor / dLiquid) %>%
             completenessCompCheck(compSetVapor) %>%
             sum())
  })
}

reMRcDVapor <- function(compSetVapor, gorVol, wtFrSystem,
                        wtFrLiquid, dLiquid, mmElemVapor, sElemVapor,
                        p = 101.325, t = 293.15, R = 8.31451) {
  return(function(dVapor) {
    checkElemNames(wtFrSystem)
    wtFrSystem %<>% checkCompStake()
    checkElemNames(wtFrLiquid)
    wtFrLiquid %<>% checkCompStake()
    compareElemNames(names(wtFrSystem), names(wtFrLiquid))
    checkElemNames(mmElemVapor)
    checkElemNames(sElemVapor)
    compareElemNames(compSetVapor, names(mmElemVapor), names(sElemVapor))
    checkRealGasP(p)
    checkRealGasT(t)

    wtFrVapor <-
      reMRcWtFrVapor(wtFrSystem, wtFrLiquid, gorVol * dVapor / dLiquid) %>%
      completenessCompCheck(compSetVapor) %>%
      abs() %>%
      normalize()

    return(wtFrVapor %>%
             wtFrToMolFr(mmElemVapor) %>%
             {
               realGasD(additiveProperty(., mmElemVapor),
                        realGasZ(., sElemVapor) %>%
                          realGasReZ(p) %T>%
                          checkRealGasZ(),
                        p, t, R)
             })
  })
}

reMRcMolFrSytem <- function(compSetVapor, gorVol, wtFrSystem, wtFrLiquid,
                            dLiquid, mmElemSystem, mmElemVapor, sElemVapor,
                            p = 101.325, t = 293.15, R = 8.31451) {
  return(function(dVapor) {
    checkElemNames(wtFrSystem)
    wtFrSystem %<>% checkCompStake()
    checkElemNames(wtFrLiquid)
    wtFrLiquid %<>% checkCompStake()
    checkElemNames(mmElemSystem)
    compareElemNames(names(wtFrSystem), names(wtFrLiquid), names(mmElemSystem))
    checkElemNames(mmElemVapor)
    checkElemNames(sElemVapor)
    compareElemNames(compSetVapor, names(mmElemVapor), names(sElemVapor))
    checkRealGasP(p)
    checkRealGasT(t)

    wtFrVapor <-
      reMRcWtFrVapor(wtFrSystem, wtFrLiquid, gorVol * dVapor / dLiquid) %>%
      completenessCompCheck(compSetVapor) %>%
      abs() %>%
      normalize()

    dVaporMRc <-
      wtFrVapor %>%
      wtFrToMolFr(mmElemVapor) %>%
      {
        realGasD(additiveProperty(., mmElemVapor),
                 realGasZ(., sElemVapor) %>%
                   realGasReZ(p) %T>%
                   checkRealGasZ(),
                 p, t, R)
      }
    gorMass <- gorVol / dLiquid * dVaporMRc
    return(mrcMolElemSystem(wtFrVapor %>%
                              completenessCompCheck(names(wtFrLiquid)),
                            wtFrLiquid,
                            gorMass,
                            mmElemSystem) %>%
             mrcMolFrSystem())
  })
}

reMRcMMSystem <- function(compSetVapor, gorVol, wtFrSystem, wtFrLiquid,
                          dLiquid, mmElemSystem, mmElemVapor, sElemVapor,
                          p = 101.325, t = 293.15, R = 8.31451) {
  return(function(dVapor) {
    checkElemNames(wtFrSystem)
    wtFrSystem %<>% checkCompStake()
    checkElemNames(wtFrLiquid)
    wtFrLiquid %<>% checkCompStake()
    checkElemNames(mmElemSystem)
    compareElemNames(names(wtFrSystem), names(wtFrLiquid), names(mmElemSystem))
    checkElemNames(mmElemVapor)
    checkElemNames(sElemVapor)
    compareElemNames(compSetVapor, names(mmElemVapor), names(sElemVapor))
    checkRealGasP(p)
    checkRealGasT(t)

    wtFrVapor <-
      reMRcWtFrVapor(wtFrSystem, wtFrLiquid, gorVol * dVapor / dLiquid) %>%
      completenessCompCheck(compSetVapor) %>%
      abs() %>%
      normalize()

    dVaporMRc <-
      wtFrVapor %>%
      wtFrToMolFr(mmElemVapor) %>%
      {
        realGasD(additiveProperty(., mmElemVapor),
                 realGasZ(., sElemVapor) %>%
                   realGasReZ(p) %T>%
                   checkRealGasZ(),
                 p, t, R)
      }
    gorMass <- gorVol / dLiquid * dVaporMRc
    return(mrcMolElemSystem(wtFrVapor %>%
                              completenessCompCheck(names(wtFrLiquid)),
                            wtFrLiquid,
                            gorMass,
                            mmElemSystem) %>%
             sum() %>%
             mrcMMbyMolSystem(gorMass))
  })
}


errorRevDiff <- function(func, valueTarget) {
  return(function(dVapor) {
    (func(dVapor) / valueTarget - 1) * 100
  })
}

errorMSE <- function(func, compTarget) {
  return(function(dVapor) {
    checkElemNames(compTarget)
    compTarget %<>% checkCompStake()
    Metrics::mse(compTarget, func(dVapor))
  })
}

derivative <- function(errorFunc, h) {
  return(function(dVapor) {
    (errorFunc(dVapor + h) - errorFunc(dVapor)) / h
  })
}

dSearch <- function(errorFunc, h = 1e-8, ...) {
  return(VGAMextra::newtonRaphson.basic(f = errorFunc,
                                        fprime = derivative(errorFunc, h),
                                        ...))
}
