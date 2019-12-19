simulation_extTADA <- function(FDR, pars0, samplenumber, rankpercentile) {
  genes = FDR$Gene
  simulated_LOF = rep(0, length(genes))
  simulated_Dmis = rep(0, length(genes))
  assignment <- sample(1:4, size = length(genes), prob = pars0[1:4,1], replace = T)
  for (i in 1:length(genes)) {
    hyperGammaMean_LGD = pars0[5,1]/(1+exp(-pars0[7,1]*(rankpercentile[i]-pars0[9,1]))) + 1;
    hyperGammaMean_Dmis = pars0[6,1]/(1+exp(-pars0[8,1]*(rankpercentile[i]-pars0[10,1]))) + 1;
    if (assignment[i] == 1) {
      simulated_LOF[i] = rpois(1, 2*samplenumber*FDR$mut_cls1[i]*rgamma(1, hyperGammaMean_LGD/pars0[11,1], pars0[11,1]))
      simulated_Dmis[i] = rpois(1, 2*samplenumber*FDR$mut_cls2[i])
    } else if (assignment[i] == 2) {
      simulated_LOF[i] = rpois(1, 2*samplenumber*FDR$mut_cls1[i])
      simulated_Dmis[i] = rpois(1, 2*samplenumber*FDR$mut_cls2[i]*rgamma(1, hyperGammaMean_Dmis/pars0[12,1], pars0[12,1]))
    } else if (assignment[i] == 3) {
      simulated_LOF[i] = rpois(1, 2*samplenumber*FDR$mut_cls1[i]*rgamma(1, hyperGammaMean_LGD/pars0[11,1], pars0[11,1]))
      simulated_Dmis[i] = rpois(1, 2*samplenumber*FDR$mut_cls2[i]*rgamma(1, hyperGammaMean_Dmis/pars0[12,1], pars0[12,1]))
    } else if (assignment[i] == 4) {
      simulated_LOF[i] = rpois(1, 2*samplenumber*FDR$mut_cls1[i])
      simulated_Dmis[i] = rpois(1, 2*samplenumber*FDR$mut_cls2[i])
    } 
  }
  result = data.frame(LOF = simulated_LOF, Dmis = simulated_Dmis)
  row.names(result) <- FDR$Gene
  return(result)
}