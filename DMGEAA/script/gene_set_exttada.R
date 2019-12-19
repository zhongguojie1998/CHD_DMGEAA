gene_set_exttada <- function(geneset, cases, samplenumber, reference, cal_pval = FALSE, outputpath=NA) {
  # form a proper input for TADA
  mytada.data = data.frame(gene.id=geneset,
                           mut.cls1=numeric(length(geneset)),
                           mut.cls2=numeric(length(geneset)),
                           mut.cls3=numeric(length(geneset)),
                           dn.cls1=numeric(length(geneset)),
                           dn.cls2=numeric(length(geneset)),
                           dn.cls3=numeric(length(geneset)))
  for (i in 1:length(geneset)) {
    index = match(as.character(geneset[i]), reference$gene_name)
    # check whether we found the gene
    if (!is.na(index)) {
      mytada.data$mut.cls1[i] = reference$p_LGD[index]
      mytada.data$mut.cls2[i] = reference$REVEL_0.5[index]
      mytada.data$mut.cls3[i] = max(reference$PHRED_20[index] - reference$REVEL_0.5[index], 0)
    } else {
      mytada.data$mut.cls1[i] = .Machine$double.eps
      mytada.data$mut.cls2[i] = .Machine$double.eps
      mytada.data$mut.cls3[i] = .Machine$double.eps
    }
    if (mytada.data$mut.cls1[i] <= 0) {
      mytada.data$mut.cls1[i] = .Machine$double.eps
    }
    if (mytada.data$mut.cls2[i] <= 0) {
      mytada.data$mut.cls2[i] = .Machine$double.eps
    }
    if (mytada.data$mut.cls3[i] <= 0) {
      mytada.data$mut.cls3[i] = .Machine$double.eps
    }
    gene_cases = cases[cases$genes==as.character(geneset[i]),]
    mytada.data$dn.cls1[i] = sum(gene_cases$classes=="LGD")
    mytada.data$dn.cls2[i] = sum(gene_cases$classes!="LGD"
                                 & gene_cases$REVEL>=0.5
                                 & !is.na(gene_cases$REVEL))
    mytada.data$dn.cls3[i] = sum(gene_cases$classes!="LGD" 
                                 & gene_cases$REVEL<0.5 
                                 & !is.na(gene_cases$REVEL) 
                                 & gene_cases$CADD>=20 
                                 & !is.na(gene_cases$CADD))
  }
  n.family = samplenumber
  n = data.frame(dn=n.family, ca=NA, cn=NA)
  sample.counts <- list(cls1=n,
                        cls2=n,
                        cls3=n)
  
  cls1.counts=data.frame(dn=mytada.data$dn.cls1, ca=NA, cn=NA)
  row.names(cls1.counts)=mytada.data$gene.id
  cls2.counts=data.frame(dn=mytada.data$dn.cls2, ca=NA, cn=NA)
  row.names(cls2.counts)=mytada.data$gene.id
  cls3.counts=data.frame(dn=mytada.data$dn.cls3, ca=NA, cn=NA)
  rownames(cls3.counts)=mytada.data$gene.id
  
  tada.counts=list(cls1=cls1.counts, 
                   cls2=cls2.counts, 
                   cls3=cls3.counts)
  mu=data.frame(cls1 = mytada.data$mut.cls1, 
                cls2 = mytada.data$mut.cls2,
                cls3 = mytada.data$mut.cls3)
  denovo.only=data.frame(cls1=TRUE,cls2=TRUE,cls3=TRUE)
  # do burden analysis to get the estimated prior
  burden_analysis = geneset_burden_analysis(geneset, cases, samplenumber, reference)
  burden1 = burden_analysis$enrichment[1]
  burden2 = burden_analysis$enrichment[4]
  burden3 = burden_analysis$enrichment[5]
  
  # use extTADA to sample values of parameters
  dataDN = mytada.data[,c('dn.cls1','dn.cls2')]
  colnames(dataDN) = c('dn_cls1','dn_cls2')
  mutRate = mytada.data[,c('mut.cls1','mut.cls2')]
  colnames(mutRate) = c('mut_cls1','mut_cls2')
  
  mouse_dev = read.table('Developmental Data/mouse_br_rnaseq3.rda.txt', header = TRUE, fill = TRUE)
  mouse_dev = mouse_dev[!is.na(mouse_dev$e14.5_rank),]
  rankPercentile = mouse_dev$e14.5_rank[match(geneset, mouse_dev$human.External.Gene.Name)]
  options(mc.cores = parallel::detectCores())
  mcmcDD <- extTADAmcmc(modelName = DNextTADA,
                        dataDN = dataDN,
                        mutRate = mutRate,
                        rankPercentile = rankPercentile,
                        Ndn = rep(samplenumber, 2),
                        nIteration = 500)

  options(repr.plot.width = 4, repr.plot.height = 3)
  par(mfrow = c(1,2))
  #plotParHeatmap(c("pi0[1]", "hyperGammaMeanDN[1]"), mcmcResult = mcmcDD)
  #plotParHeatmap(c("pi0[2]", "hyperGammaMeanDN[2]"), mcmcResult = mcmcDD)
  
  pars0 = estimatePars(pars = mcmcDD@sim$fnames_oi,
                       mcmcResult = mcmcDD)
  parsFDR <- list(percentileA = pars0[,1][5:6],
                  percentileB = pars0[,1][7:8],
                  percentileC = pars0[,1][9:10],
                  pi0 = pars0[,1][1:4],
                  nfamily = rep(samplenumber, 2),
                  betaPars = c(6.7771073, -1.7950864, -0.2168248),
                  rankPercentile = rankPercentile)

  dataFDR <- calculateFDR(pars = parsFDR,
                          dnData = dataDN,
                          mutData = mutRate,
                          geneName = geneset)
  result = list("mcmcDD"=mcmcDD, "pars0"=pars0, "dataFDR"=dataFDR)
  return(result)
}
