source("geneset_burden_analysis.R")
library(denovolyzeR)
library(readxl)
mydata = read_excel("PCGC_DNVs_Freeze201907.xlsx", sheet = "DNVs", na = ".")
trois = read_excel("PCGC_DNVs_Freeze201907.xlsx", sheet = "Trios", na = ".")
samplenumber = length(unique(trois$IID))
effs = mydata$GeneEff
standards = c("frameshift"="LGD",
              "stop_gained"="LGD",
              "missense"="mis",
              "splice_region"=NA,
              "synonymous"="syn",
              "splice_acceptor"="LGD",
              "splice_donor"="LGD",
              "start_lost"="LGD",
              "inframe_deletion"=NA,
              "inframe_insertion"=NA,
              "synonymous;missense"="mis",
              "stop_lost"="LGD",
              "protein_altering"="LGD",
              "stop_retained"="syn")
standard_var = effs
topick = c()
for (i in 1:length(effs)) {
  standard_var[i] = standards[effs[i]]
  if (!is.na(standards[effs[i]])) {
    topick = c(topick, i)
  }
}
## first do all genes
reference = read.table("mutationrate.20190816.canonical.map0.txt", sep = "\t", header = TRUE)

#here I use genes in the ExAC data set as all genes
ExAC = read_excel("ExAC/nature19057-SI Table 13.xlsx", sheet = "Gene Constraint")
cases = data.frame(genes = mydata$Symbol[topick], classes = standard_var[topick], REVEL = mydata$REVEL[topick], CADD13 = mydata$CADD13[topick])

# only use genes in the reference dataset
geneset = unique(reference$gene_name)
alldev = read.table('Developmental Data/mouse_br_rnaseq3.rda.txt', header = TRUE, fill = TRUE)
alldev = alldev[!is.na(alldev$e14.5_rank)&!is.na(alldev$human.External.Gene.Name),]
alldev_geneset = alldev$human.External.Gene.Name
geneset = unique(alldev_geneset[alldev_geneset %in% geneset])

## 4 pi0s + gamma function model
DMGEAA = dir('DMGEAA/script/', '.R$')
for (ii in DMGEAA) {
  source(paste0('DMGEAA/script/', ii))
}
# only use genes in the developmental gene dataset
alldev_geneset = geneset
# alldev_geneset_exttada includes the mcmc samples, FDR of genes, and the estimated pars
alldev_geneset_exttada = gene_set_exttada(alldev_geneset, cases, samplenumber, reference)
rank_percentile = alldev$e14.5_rank[match(alldev_geneset, alldev$human.External.Gene.Name)]
# compare with simulation
alldev_simulation_compare = compare_with_simulation(alldev_geneset_exttada, samplenumber = samplenumber, rank_percentile = rank_percentile)

# alldev_minus deleted three well known genes
alldev_geneset_minus = alldev_geneset[-which(alldev_geneset=="KMT2D" | alldev_geneset=="CHD7" | alldev_geneset=="PTPN11")]
# alldev_geneset_minus_exttada includes the mcmc samples, FDR of genes, and the estimated pars
alldev_geneset_minus_exttada = gene_set_exttada(alldev_geneset_minus, cases, samplenumber, reference)
rank_percentile = alldev$e14.5_rank[match(alldev_geneset_minus, alldev$human.External.Gene.Name)]
# compare with simulation
alldev_minus_simulation_compare = compare_with_simulation(alldev_geneset_minus_exttada, samplenumber = samplenumber, rank_percentile = rank_percentile)


