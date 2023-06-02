#!/usr/bin/env Rscript

library(hyprcoloc)
library(argparse)
library(data.table)
library(IGUtilityPackage)

parser <- ArgumentParser(description = 'Run HyprColoc for every locus to detect colocalisation between pQTL and eQTL datasets.')
parser$add_argument('--harmonised_data', metavar = 'file', type = 'character',
                    help = 'R data file which contains list with harmonised sumstats matrices')
parser$add_argument('--gene_id', type = 'character',
                    help = 'Gene name for which the colocalisation is done.')
parser$add_argument('--protein_id', type = 'character',
                    help = 'Protein ID for which the colocalisation is done.')

args <- parser$parse_args()

harmonised_data <- readRDS(args$harmonised_data)

eqtl <- harmonised_data$sumstats1
pqtl <- harmonised_data$sumstats2

res <- data.table(iteration = NA, 
traits = NA, 
posterior_prob = NA, 
regional_prob = NA, 
candidate_snp = NA, 
posterior_explained_by_snp = NA, 
dropped_trait = NA,
nr_snps_included = NA,
analysis = NA, 
lead_variant = NA,
region = NA,
dataset1_effect = NA,
dataset2_effect = NA,
dataset2_p = NA,
dataset2_min_p = NA)
res <- res[-1, ]

#######################
# Analyse pQTL regions#
#######################

# Find lead variants for pQTLs
message("Analyse pQTL-eQTL")
pqtl_lead <- IdentifyLeadSNPs(data = pqtl,
snp_id_col = "SnpId",
snp_chr_col = "chr",
snp_pos_col = "pos",
eff_all_col = "ea",
other_all_col = "nea",
beta_col = "beta",
se_col = "se")

for (i in 1:nrow(pqtl_lead)){
pqtl_reg <- pqtl[chr == pqtl_lead$chr[i] & pos > pqtl_lead$pos[i] - 1000000 & pos < pqtl_lead$pos[i] + 1000000]
eqtl_reg <- eqtl[SnpId %in% pqtl_reg$SnpId]

res_temp <- hyprcoloc(as.matrix(data.table(pqtl_reg$beta, eqtl_reg$beta)), as.matrix(data.table(pqtl_reg$se, eqtl_reg$se)), 
trait.names = c("pQTL", "eQTL"))  

res_temp <- res_temp$results
res_temp$nr_snps_included <- nrow(pqtl_reg)

res_temp$analysis = "pQTL-eQTL"
res_temp$lead_variant <- pqtl_lead$SNP[i]
res_temp$region <- paste0(pqtl_lead$chr[i], ":", pqtl_lead$pos[i] - 1000000, "-", pqtl_lead$pos[i] + 1000000)

res_temp$dataset1_effect = pqtl_lead$beta[i]
res_temp$dataset2_effect = eqtl_reg[eqtl_reg$SnpId == pqtl_lead$SNP[i], ]$beta

eqtl_reg$P <- ZtoP(eqtl_reg$beta/eqtl_reg$se)

res_temp$dataset2_p <- eqtl_reg[eqtl_reg$SnpId == pqtl_lead$SNP[i], ]$P
res_temp$dataset2_min_p <- min(eqtl_reg$P)

if (!is.na(res_temp$candidate_snp)){res_temp$candidate_snp <- eqtl_reg$SnpId[as.numeric(res_temp$candidate_snp)]}


message(paste0("Analysed locus ", i, "/", nrow(pqtl_lead)))

res <- rbind(res, res_temp)
}

#######################
# Analyse eQTL regions#
#######################
message("Analyse eQTL-pQTL")
# Find lead variants for eQTLs
eqtl_lead <- IdentifyLeadSNPs(data = eqtl,
snp_id_col = "SnpId",
snp_chr_col = "chr",
snp_pos_col = "pos",
eff_all_col = "ea",
other_all_col = "nea",
beta_col = "beta",
se_col = "se")

for (i in 1:nrow(eqtl_lead)){
pqtl_reg <- pqtl[chr == eqtl_lead$chr[i] & pos > eqtl_lead$pos[i] - 1000000 & pos < eqtl_lead$pos[i] + 1000000]
eqtl_reg <- eqtl[SnpId %in% pqtl_reg$SnpId]

res_temp <- hyprcoloc(as.matrix(data.table(pqtl_reg$beta, eqtl_reg$beta)), as.matrix(data.table(pqtl_reg$se, eqtl_reg$se)), 
trait.names = c("pQTL", "eQTL"))  

res_temp <- res_temp$results
res_temp$nr_snps_included <- nrow(pqtl_reg)

res_temp$analysis = "eQTL-pQTL"
res_temp$lead_variant <- eqtl_lead$SNP[i]
res_temp$region <- paste0(pqtl_lead$chr[i], ":", pqtl_lead$pos[i] - 1000000, "-", pqtl_lead$pos[i] + 1000000)

res_temp$dataset1_effect = eqtl_lead$beta[i]
res_temp$dataset2_effect = pqtl_reg[pqtl_reg$SnpId == eqtl_lead$SNP[i], ]$beta

pqtl_reg$P <- ZtoP(pqtl_reg$beta/pqtl_reg$se)

res_temp$dataset2_p <- pqtl_reg[pqtl_reg$SnpId == eqtl_lead$SNP[i], ]$P
res_temp$dataset2_min_p <- min(pqtl_reg$P)

if (!is.na(res_temp$candidate_snp)){res_temp$candidate_snp <- eqtl_reg$SnpId[as.numeric(res_temp$candidate_snp)]}

message(paste0("Analysed locus ", i, "/", nrow(eqtl_lead)))

res <- rbind(res, res_temp)
}

res <- data.table(gene = args$gene_id, protein = args$protein_id, res)

fwrite(res, paste0(args$protein_id, "__", args$gene_id, ".txt"), sep = "\t")

