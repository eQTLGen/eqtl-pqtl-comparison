#!/usr/bin/env Rscript

library(TwoSampleMR)
library(data.table)
library(IGUtilityPackage)
library(argparse)

setDTthreads(1)

parser <- ArgumentParser(description = 'Run Genome-wide two-sample MR between two datasets.')

parser$add_argument('--harmonised_data', metavar = 'file', type = 'character',
                    help = 'R data file which contains list with harmonised sumstats matrices.')
parser$add_argument('--gene_id', type = 'character',
                    help = 'Gene name for which the MR is done.')
parser$add_argument('--protein_id', type = 'character',
                    help = 'Protein ID for which the MR is done.')

args <- parser$parse_args()


and <- readRDS(args$harmonised_data)

eqtl <- and$sumstats1
pqtl <- and$sumstats2

# gene -> protein
message("Analysis: gene -> protein")  

eqtl_lead <- IdentifyLeadSNPs(
       eqtl,
       window = 1e+06,
       Pthresh = 5e-08,
       snp_id_col = "SnpId",
       snp_chr_col = "chr",
       snp_pos_col = "pos",
       eff_all_col = "ea",
       other_all_col = "nea",
       beta_col = "beta",
       se_col = "se",
       p_col = NULL
     )

message("gene -> protein: lead variants found")
   
exposure <- eqtl_lead[, c(1, 6, 7, 4, 5)]
colnames(exposure) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure")
exposure$`eaf.exposure` <- 0.2
exposure$exposure <- args$gene_id
exposure$`id.exposure` <- args$gene_id


outcome <- pqtl[, c(1, 6, 7, 4, 5)]
colnames(outcome) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome")
outcome$`eaf.outcome` <- 0.2
outcome$outcome <- args$protein_id
outcome$`id.outcome` <- args$protein_id

inp <- harmonise_data(exposure, outcome, action = 1)
message("gene -> protein: data harmonised")

res_eqtl_pqtl <- mr(inp)
message("gene -> protein: MR finished")

# protein -> gene
message("Analysis: protein -> gene")  
pqtl_lead <- IdentifyLeadSNPs(
       pqtl,
       window = 1e+06,
       Pthresh = 5e-08,
       snp_id_col = "SnpId",
       snp_chr_col = "chr",
       snp_pos_col = "pos",
       eff_all_col = "ea",
       other_all_col = "nea",
       beta_col = "beta",
       se_col = "se",
       p_col = NULL
     )

message("protein -> gene: lead variants found")
exposure <- pqtl_lead[, c(1, 6, 7, 4, 5)]
colnames(exposure) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure")
exposure$`eaf.exposure` <- 0.2
exposure$exposure <- args$protein_id
exposure$`id.exposure` <- args$protein_id

outcome <- eqtl[, c(1, 6, 7, 4, 5)]
colnames(outcome) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome")
outcome$`eaf.outcome` <- 0.2
outcome$outcome <- args$gene_id
outcome$`id.outcome` <- args$gene_id

inp <- harmonise_data(exposure, outcome, action = 1)
message("protein -> gene: data harmonised")
res_pqtl_eqtl <- mr(inp)
message("protein -> gene: MR finished")

head(res_eqtl_pqtl)
head(res_pqtl_eqtl)

res <- rbind(res_eqtl_pqtl, res_pqtl_eqtl)

fwrite(res, paste0(args$protein_id, "_", args$gene_id, ".txt"), sep = "\t", quote = FALSE)
message("Analysis finished!")

