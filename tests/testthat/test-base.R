testthat::test_that("normalize test", {
  set.seed(1)
  testthat::expect_equal(sum(normalize(runif(10, min = 0, max = 1))), 1)
  testthat::expect_equal(sum(normalize(runif(10, min = 0, max = 100))), 1)
})

testthat::test_that("molFrToWtFr test", {
  molFr <- c(0.071, 0.785, 0.045, 0.036, 0.016, 0.031, 0.004, 0.004, 0.008) %>% normalize()
  names(molFr) <- c("N2", "C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6")
  MMelem <- c(28, 16, 30.1, 44.1, 58.1, 58.1, 72.1, 72.1, 86.2)
  names(MMelem) <- c("N2", "C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6")
  wtFrTarget <- c(0.0925, 0.5846, 0.0630, 0.0739, 0.0433, 0.0838, 0.0134, 0.0134, 0.0321)
  names(wtFrTarget) <- c("N2", "C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6")
  testthat::expect_equal(molFrToWtFr(molFr, MMelem) %>% round(2), wtFrTarget %>% round(2))
})

testthat::test_that("wtFrToMolFr test", {
  wtFr <- c(0.0925, 0.5846, 0.0630, 0.0739, 0.0433, 0.0838, 0.0134, 0.0134, 0.0321) %>% normalize()
  names(wtFr) <- c("N2", "C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6")
  MMelem <- c(28, 16, 30.1, 44.1, 58.1, 58.1, 72.1, 72.1, 86.2)
  names(MMelem) <- c("N2", "C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6")
  molFrTarget <- c(0.071, 0.785, 0.045, 0.036, 0.016, 0.031, 0.004, 0.004, 0.008)
  names(molFrTarget) <- c("N2", "C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6")
  testthat::expect_equal(wtFrToMolFr(wtFr, MMelem) %>% round(3), molFrTarget)
})

testthat::test_that("additiveProperty test", {
  molFr <- c(0.071, 0.785, 0.045, 0.036, 0.016, 0.031, 0.004, 0.004, 0.008) %>% normalize()
  names(molFr) <- c("N2", "C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6")
  MMelem <- c(28, 16, 30.1, 44.1, 58.1, 58.1, 72.1, 72.1, 86.2)
  names(MMelem) <- c("N2", "C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6")
  MMtarget <- 21.49
  testthat::expect_equal(additiveProperty(molFr, MMelem) %>% round(2), MMtarget)
})

testthat::test_that("idelGasVM test", {
  p <- 101.325
  t <- 293.15
  R <- 8.31451
  GasVMtarget <- 24.05525
  testthat::expect_equal(idealGasVM(p, t, R) %>% round(5), GasVMtarget)
})

testthat::test_that("idelGasG test", {
  MMsubstance <- 28.69570
  MMreference <- 28.96260
  GasGtarget <- 0.9908
  testthat::expect_equal(idealGasG(MMsubstance, MMreference) %>% round(4), GasGtarget)
})

testthat::test_that("idelGasD test", {
  MM <- 28.69570
  p <- 101.325
  t <- 293.15
  R <- 8.31451
  GasDtarget <- 1.1929
  testthat::expect_equal(idealGasD(MM, p, t, R) %>% round(4), GasDtarget)
})

testthat::test_that("realGasZ test", {
  molFr <- c(0.933212, 0.025656, 0.015368, 0.010350, 0.015414)
  names(molFr) <- c("C1", "C2", "C3", "N2", "CO2")
  sComponent <- c(0.04452, 0.09190, 0.13440, 0.01700, 0.07520)
  names(sComponent) <- c("C1", "C2", "C3", "N2", "CO2")
  Ztarget <- 0.99776224
  testthat::expect_equal(realGasZ(molFr, sComponent) %>% round(4), Ztarget %>% round(4))
})

testthat::test_that("realGasReZ test", {
  Z <- 0.9930
  p <- 110
  Ztarget <- 0.9924
  testthat::expect_equal(realGasReZ(Z, p) %>% round(4), Ztarget)
})

testthat::test_that("realGasVM test", {
  Z <- 0.99776224
  p <- 101.325
  t <- 288.15
  R <- 8.314462
  VMtarget <- 23.591917
  testthat::expect_equal(realGasVM(Z, p, t, R) %>% round(5), VMtarget %>% round(5))
})

testthat::test_that("realGasG test", {
  MMsubstance <- 28.69570
  MMreference <- 28.96260
  Zsubstance <- 0.9930
  Zreference <- 0.99963
  Gtarget <- 0.9974
  testthat::expect_equal(realGasG(MMsubstance, Zsubstance, MMreference, Zreference) %>% round(4), Gtarget %>% round(4))
})

testthat::test_that("realGasD test", {
  MM <- 28.69570
  Z <- 0.9930
  p <- 101.325
  t <- 293.15
  R <- 8.314462
  Dtarget <- 1.2013
  testthat::expect_equal(realGasD(MM, Z, p, t, R) %>% round(4), Dtarget)
})

testthat::test_that("mrcMolElemSystem test", {
  wtFrVapor <- c(0.09280, 0.388760, 0.14681, 0.25030, 0.07900, 0.02510, 0.00690, 0.01150) %>% normalize()
  names(wtFrVapor) <- c("N2", "C1", "C2", "C3", "C4", "C5", "C6", "C7+")
  wtFrLiquid <- c(0, 0.00021, 0.00132, 0.01176, 0.01423, 0.01816, 0.02268, 0.93164) %>% normalize()
  names(wtFrLiquid) <- c("N2", "C1", "C2", "C3", "C4", "C5", "C6", "C7+")
  gorMass <- 0.2236
  MMelemSystem <- c(28, 16, 30.1, 44.1, 58.1, 72.1, 86.2, 261.8)
  names(MMelemSystem) <- c("N2", "C1", "C2", "C3", "C4", "C5", "C6", "C7+")
  molSystemTarget <- c(0.00074, 0.00543, 0.00113, 0.00154, 0.00055, 0.00033, 0.00028, 0.00357)
  names(molSystemTarget) <- c("N2", "C1", "C2", "C3", "C4", "C5", "C6", "C7+")
  testthat::expect_equal(mrcMolElemSystem(wtFrVapor, wtFrLiquid, gorMass, MMelemSystem) %>% round(3), molSystemTarget %>% round(3))
})

testthat::test_that("mrcWtFrSystem test", {
  wtFrVapor <- c(0.09280, 0.388760, 0.14681, 0.25030, 0.07900, 0.02510, 0.00690, 0.01150) %>% normalize()
  names(wtFrVapor) <- c("N2", "C1", "C2", "C3", "C4", "C5", "C6", "C7+")
  wtFrLiquid <- c(0, 0.00021, 0.00132, 0.01176, 0.01423, 0.01816, 0.02268, 0.93164) %>% normalize()
  names(wtFrLiquid) <- c("N2", "C1", "C2", "C3", "C4", "C5", "C6", "C7+")
  gorMass <- 0.2236
  wtFrSystemTarget <- c(0.01696, 0.07101, 0.02790, 0.05535, 0.02606, 0.01943, 0.01979, 0.76350)
  names(wtFrSystemTarget) <- c("N2", "C1", "C2", "C3", "C4", "C5", "C6", "C7+")
  testthat::expect_equal(mrcWtFrSystem(wtFrVapor, wtFrLiquid, gorMass) %>% round(3), wtFrSystemTarget %>% round(3))
})

testthat::test_that("mrcMolFrSystem test", {
  molSystem <- c(0.00074, 0.00543, 0.00113, 0.00154, 0.00055, 0.00033, 0.00028, 0.00357)
  names(molSystem) <- c("N2", "C1", "C2", "C3", "C4", "C5", "C6", "C7+")
  molFrSystmeTarget <- c(0.05453, 0.40072, 0.08347, 0.11302, 0.04041, 0.02428, 0.02070, 0.26287)
  names(molFrSystmeTarget) <- c("N2", "C1", "C2", "C3", "C4", "C5", "C6", "C7+")
  testthat::expect_equal(mrcMolFrSystem(molSystem) %>% round(2), molFrSystmeTarget %>% round(2))
})

testthat::test_that("mrcMMbyMolSystem test", {
  gorMass <- 0.2236
  molSystem <- 0.01357
  mmSystemTarget <- 90.2
  testthat::expect_equal(mrcMMbyMolSystem(molSystem, gorMass) %>% round(1), mmSystemTarget)
})

testthat::test_that("reMRcWtFrVapor test", {
  compSet <- c("CO2", "C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6", "C7",
               "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16",
               "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25",
               "C26", "C27", "C28", "C29", "C30", "C31", "C32", "C33", "C34",
               "C35", "C36+")
  compSetVapor <- c("CO2", "C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5",
                    "C6", "C7", "C8", "C9")
  wtFrSystem <- c(0.0017, 0.9254, 0.7209, 1.4776, 1.0429, 1.7663, 1.4162,
                  1.3952, 1.9360, 3.9130, 6.2578, 3.7748, 3.5375, 2.9115,
                  2.7836, 3.0401, 3.0256, 3.5393, 3.5174, 4.8179, 3.9671,
                  3.0089, 3.5876, 3.4420, 3.1631, 2.9597, 2.7074, 2.7246,
                  2.2701, 2.2280, 1.7873, 1.7732, 1.4928, 1.1937, 0.9615,
                  0.8617, 0.7310, 0.5845, 8.7551) %>% normalize()
  names(wtFrSystem) <- compSet
  wtFrLiquid <- c(0, 0, 0.0592, 0.3767, 0.6209, 1.2669, 1.2574, 1.2973, 1.8535,
                  3.9684, 6.5006, 3.9460, 3.7008, 3.0459, 2.9122, 3.1805, 3.1654,
                  3.7028, 3.6799, 5.0404, 4.1503, 3.1478, 3.7533, 3.6009, 3.3091,
                  3.0964, 2.8324, 2.8504, 2.3749, 2.3309, 1.8699, 1.8550, 1.5617,
                  1.2488, 1.0059, 0.9015, 0.7648, 0.6115, 9.1594) %>% normalize()
  names(wtFrLiquid) <- compSet
  dLiquid <- 817.531
  gorVol <- 26.723622049749
  dVapor <- 1.41277784469918
  wtFrVaporTarget <- c(0.0378, 20.9628, 15.0493, 25.3160, 10.1791, 12.5809,
                       4.8553, 3.5159, 3.7226, 2.7135, 1.0007, 0.0661) %>% normalize()
  names(wtFrVaporTarget) <- compSetVapor
  testthat::expect_equal(reMRcWtFrVapor(wtFrSystem, wtFrLiquid, gorVol * dVapor / dLiquid)[1:12] %>% round(3), wtFrVaporTarget %>% round(3))
})

testthat::test_that("dVaporMB test", {
  VF <- 1.06973923097113
  gorVol <- 26.723622049749
  dSystem <- 800.737
  dLiquid <- 817.531
  dVaporTarget <- 1.4612
  testthat::expect_equal(dVaporMB(VF, gorVol, dSystem, dLiquid) %>% round(4), dVaporTarget)
})
