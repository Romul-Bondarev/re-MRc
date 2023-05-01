testthat::test_that("checkElemNames test", {
  xFr <- c(0.50, 0.25, 0, 0, 0, 0.25, 0)
  testthat::expect_error(checkElemNames(xFr))
})

testthat::test_that("compareElemNames test", {
  xFrNames <- c("C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5")
  testthat::expect_error(compareElemNames(xFrNames))

  yFrNames <- c("N2", "CO2", "C1", "C2", "C3", "IC4", "NC4")
  testthat::expect_error(compareElemNames(xFrNames, yFrNames))
})

testthat::test_that("completenessCompCheck test", {
  xFr <- c(0.50, 0.25, 0, 0, 0, 0.25, 0)
  names(xFr) <- c("C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5")

  compNames <- c("C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5", "C6", "C7", "C8", "C9")
  testthat::expect_equal(names(completenessCompCheck(xFr, compNames)), compNames)

  compNames <- c("C1", "C2", "IC5", "C6", "C7", "C8", "C9")
  testthat::expect_equal(names(completenessCompCheck(xFr, compNames)), compNames)

  yFr <- c(0, 0.25, 0, 0, 0, 0.25, 0.50)
  names(yFr) <- c("NC5", "IC5", "NC4", "IC4", "C3", "C2", "C1")
  testthat::expect_equal(names(completenessCompCheck(yFr, compNames)), compNames)
})

testthat::test_that("checkCompStake test", {
  xFr <- c(55, 25, 25, 0, 0, 0, -5)
  names(xFr) <- c("C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5")
  testthat::expect_error(checkCompStake(xFr))

  xFr <- c(55, 25, 25, 0, 0, 0, 0)
  names(xFr) <- c("C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5")
  testthat::expect_error(checkCompStake(xFr))

  xFr <- c(50, 25, 0, 0, 0, 25, 0)
  names(xFr) <- c("C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5")
  testthat::expect_message(testthat::expect_equal(sum(checkCompStake(xFr)), 1))

  xFr <- c(0.50, 0.25, 0, 0, 0, 0.25, 0)
  names(xFr) <- c("C1", "C2", "C3", "IC4", "NC4", "IC5", "NC5")
  testthat::expect_identical(sum(checkCompStake(xFr)), 1)
})

testthat::test_that("checkRealGasZ, checkRealGasP, checkRealGasT test", {
  testthat::expect_warning(checkRealGasZ(0.9))
  testthat::expect_warning(checkRealGasZ(1.1))

  testthat::expect_warning(checkRealGasP(85))
  testthat::expect_warning(checkRealGasP(115))

  testthat::expect_warning(checkRealGasT(270))
})

testthat::test_that("reMRcSumVapor test", {
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
  mmElemSystem <- c(44.01, 16.043, 30.07, 44.097, 58.123, 58.123, 72.15, 72.15,
                    84, 96, 107, 121, 134, 147, 161, 175, 190, 206, 222, 237, 251,
                    263, 275, 291, 305, 318, 331, 345, 359, 374, 388, 402, 416,
                    430, 444, 458, 472, 486, 608.35)
  names(mmElemSystem) <- compSet
  sElemVapor <- c(0.0728, 0.0436, 0.0894, 0.1288, 0.1703, 0.1783, 0.2168, 0.2345,
                  0.2846, 0.3521, 0.4278, 0.5148)
  names(sElemVapor) <- compSetVapor
  dLiquid <- 817.531
  gorVol <- 26.723622049749
  dVapor <- 1.41277784469918

  sumVapor <- reMRcSumVapor(compSetVapor, gorVol, wtFrSystem, wtFrLiquid, dLiquid)
  testthat::expect_equal(sumVapor(dVapor) %>% round(4), 1 %>% round(4))
})

testthat::test_that("reMRcDVapor test", {
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
  mmElemSystem <- c(44.01, 16.043, 30.07, 44.097, 58.123, 58.123, 72.15, 72.15,
                    84, 96, 107, 121, 134, 147, 161, 175, 190, 206, 222, 237, 251,
                    263, 275, 291, 305, 318, 331, 345, 359, 374, 388, 402, 416,
                    430, 444, 458, 472, 486, 608.35)
  names(mmElemSystem) <- compSet
  sElemVapor <- c(0.0728, 0.0436, 0.0894, 0.1288, 0.1703, 0.1783, 0.2168, 0.2345,
                  0.2846, 0.3521, 0.4278, 0.5148)
  names(sElemVapor) <- compSetVapor
  dLiquid <- 817.531
  gorVol <- 26.723622049749
  dVapor <- 1.41277784469918

  dVaporCurrent <- reMRcDVapor(compSetVapor, gorVol, wtFrSystem, wtFrLiquid, dLiquid, mmElemSystem %>% completenessCompCheck(compSetVapor), sElemVapor)
  testthat::expect_equal(dVaporCurrent(dVapor) %>% round(4), dVapor %>% round(4))
})

testthat::test_that("reMRcMolFrSytem test", {
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
  molFrSystemTarget <- c(0.0059, 8.9739, 3.7290, 5.2120, 2.7911, 4.7274, 3.0533,
                   3.0081, 3.5851, 6.3404, 9.0974, 4.8527, 4.1064, 3.0809,
                   2.6894, 2.7023, 2.4771, 2.6726, 2.4646, 3.1622, 2.4585,
                   1.7796, 2.0293, 1.8399, 1.6132, 1.4478, 1.2723, 1.2285,
                   0.9836, 0.9267, 0.7166, 0.6861, 0.5582, 0.4318, 0.3369,
                   0.2927, 0.2409, 0.1871, 2.2386) %>% normalize()
  names(molFrSystemTarget) <- compSet
  wtFrLiquid <- c(0, 0, 0.0592, 0.3767, 0.6209, 1.2669, 1.2574, 1.2973, 1.8535,
                  3.9684, 6.5006, 3.9460, 3.7008, 3.0459, 2.9122, 3.1805, 3.1654,
                  3.7028, 3.6799, 5.0404, 4.1503, 3.1478, 3.7533, 3.6009, 3.3091,
                  3.0964, 2.8324, 2.8504, 2.3749, 2.3309, 1.8699, 1.8550, 1.5617,
                  1.2488, 1.0059, 0.9015, 0.7648, 0.6115, 9.1594) %>% normalize()
  names(wtFrLiquid) <- compSet
  mmElemSystem <- c(44.01, 16.043, 30.07, 44.097, 58.123, 58.123, 72.15, 72.15,
                    84, 96, 107, 121, 134, 147, 161, 175, 190, 206, 222, 237, 251,
                    263, 275, 291, 305, 318, 331, 345, 359, 374, 388, 402, 416,
                    430, 444, 458, 472, 486, 608.35)
  names(mmElemSystem) <- compSet
  sElemVapor <- c(0.0728, 0.0436, 0.0894, 0.1288, 0.1703, 0.1783, 0.2168, 0.2345,
                  0.2846, 0.3521, 0.4278, 0.5148)
  names(sElemVapor) <- compSetVapor
  dLiquid <- 817.531
  gorVol <- 26.723622049749
  dVapor <- 1.41277784469918

  molFrSystem <- reMRcMolFrSytem(compSetVapor, gorVol, wtFrSystem, wtFrLiquid, dLiquid, mmElemSystem, mmElemSystem %>% completenessCompCheck(compSetVapor), sElemVapor)
  testthat::expect_equal(molFrSystem(dVapor) %>% round(4), molFrSystemTarget %>% round(4))
})
