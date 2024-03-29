library(denovolyzeR)
library(readxl)
# read trio and mutation data
mydata = read_excel("PCGC_DNVs_Freeze201907.xlsx", sheet = "DNVs", na = ".")
trois = read_excel("PCGC_DNVs_Freeze201907.xlsx", sheet = "Trios", na = ".")
samplenumber = length(unique(trois$IID))
effs = mydata$GeneEff
# do mutation annotation
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

# here I use genes in the ExAC data set as all genes
ExAC = read_excel("ExAC/nature19057-SI Table 13.xlsx", sheet = "Gene Constraint")
# get mutation cases
cases = data.frame(genes = mydata$Symbol[topick], classes = standard_var[topick], REVEL = mydata$REVEL[topick], CADD13 = mydata$CADD13[topick])

# only use genes in the reference dataset
geneset = unique(reference$gene_name)
# gene expression data in developmental process
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
# alldev_geneset_DMGEAA includes the mcmc samples, FDR of genes, and the estimated pars
# for pi, pi0[1],pi0[2],pi0[3],pi0[4] correspond to \pi_1, \pi_2, \pi_3 and \pi_0 in the paper
# for percentileA, percentileB, percentileC correspond to A, B, C in the paper, and [1] is LGD, [2] is Dmis
alldev_geneset_DMGEAA = gene_set_DMGEAA(alldev_geneset, cases, samplenumber, reference)
rank_percentile = alldev$e14.5_rank[match(alldev_geneset, alldev$human.External.Gene.Name)]
# compare with simulation
alldev_simulation_compare = compare_with_simulation(alldev_geneset_DMGEAA, samplenumber = samplenumber, rank_percentile = rank_percentile)

# alldev_minus deleted three well known genes
alldev_geneset_minus = alldev_geneset[-which(alldev_geneset=="KMT2D" | alldev_geneset=="CHD7" | alldev_geneset=="PTPN11")]
# alldev_geneset_minus_DMGEAA includes the mcmc samples, FDR of genes, and the estimated pars
# for pi, pi0[1],pi0[2],pi0[3],pi0[4] correspond to \pi_1, \pi_2, \pi_3 and \pi_0 in the paper
# for percentileA, percentileB, percentileC correspond to A, B, C in the paper, and [1] is LGD, [2] is Dmis
alldev_geneset_minus_DMGEAA = gene_set_DMGEAA(alldev_geneset_minus, cases, samplenumber, reference)
rank_percentile = alldev$e14.5_rank[match(alldev_geneset_minus, alldev$human.External.Gene.Name)]
# compare with simulation
alldev_minus_simulation_compare = compare_with_simulation(alldev_geneset_minus_DMGEAA, samplenumber = samplenumber, rank_percentile = rank_percentile)


