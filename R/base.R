normalize <- function(xFr) {
  return(xFr / sum(xFr))
}

molFrToWtFr <- function(molFr, mmElem) {
  return((molFr * mmElem) %>% normalize())
}

wtFrToMolFr <- function(wtFr, mmElem) {
  return((wtFr / mmElem) %>% normalize())
}

additiveProperty <- function(xFr, propElem) {
  xFr %<>% matrix(nrow = 1)
  propElem %<>% matrix(ncol = 1)
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

realGasVM <- function(z, p, t, R) {
  return(z * idealGasVM(p, t, R))
}

realGasG <- function(mmSubstance, zSubstance, mmReference, zReference) {
  return(idealGasG(mmSubstance, mmReference) * (zReference / zSubstance))
}

realGasD <- function(mm, z, p, t, R) {
  return(idealGasD(mm, p, t, R) / z)
}

mrcMolElemSystem <- function(wtFrVapor, wtFrLiquid, gorMass, mmElemSystem) {
  return((wtFrVapor * gorMass + wtFrLiquid) / mmElemSystem)
}

mrcWtFrSystem <- function(wtFrVapor, wtFrLiquid, gorMass) {
  return((wtFrVapor * gorMass + wtFrLiquid) / (1 + gorMass))
}

mrcMolFrSystem <- function(molElemSystem, molSystem=sum(molElemSystem)) {
  return(molElemSystem / molSystem)
}

mrcMMbyMolSystem <- function(molSystem, gorMass) {
  return((1 + gorMass) / molSystem)
}

reMRcWtFrVapor <- function(wtFrSystem, wtFrLiquid, gorMass) {
  return(wtFrSystem + (wtFrSystem - wtFrLiquid) / gorMass)
}

dVaporMB <- function(VF, gorVol, dSystem, dLiquid) {
  return((VF * dSystem - dLiquid) / gorVol)
}
