compare_with_simulation <- function(geneset_exttada, samplenumber, rank_percentile = NULL) {
  # do simulations
  observation = geneset_exttada$dataFDR[c("dn_cls1", "dn_cls2")]
  row.names(observation) <- geneset_exttada$dataFDR$Gene
  LGD_compare = matrix(0, 26, 2)
  Dmis_compare = matrix(0, 26, 2)
  for (i in 0:25) {
    LGD_compare[i+1, 1] = sum(observation[,1]==i)
    Dmis_compare[i+1, 1] = sum(observation[,2]==i)
  }
  simulationN = 100
  for (j in 1:simulationN) {
    simulation = simulation_extTADA(geneset_exttada$dataFDR,geneset_exttada$pars0,samplenumber, rankpercentile = rank_percentile)
    for (i in 0:25) {
      LGD_compare[i+1, 2] = LGD_compare[i+1, 2]+sum(simulation[,1]==i)
      Dmis_compare[i+1, 2] = Dmis_compare[i+1, 2]+sum(simulation[,2]==i)
    }
  }
  LGD_compare[, 2]=LGD_compare[, 2]/simulationN
  Dmis_compare[, 2]=Dmis_compare[, 2]/simulationN
  fisher.test(LGD_compare, workspace = 2e8)
  fisher.test(Dmis_compare, workspace = 2e8)
  result = list(LGD_compare = LGD_compare,
                Dmis_compare = Dmis_compare)
}