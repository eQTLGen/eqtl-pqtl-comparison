#!/usr/bin/env Rscript

# Load libraries
library(IGUtilityPackage)
library(argparse)
library(data.table)
library(stringr)

setDTthreads(1)

parser <- ArgumentParser(description = 'Harmonise eQTL and pQTL datasets.')
parser$add_argument('--pqtl_folder', type = 'character',
                    help = 'A tab-delimited file with pQTL results.')
parser$add_argument('--eqtl_folder', type = 'character',
                    help = 'A tab-delimited file with eQTL results.')
parser$add_argument('--allele_info', metavar = 'file', type = 'character',
                    help = 'A tab-delimited file with allele information.')
parser$add_argument('--gene_id', type = 'character',
                    help = 'Gene name for which the colocalisation is done.')
parser$add_argument('--protein_id', type = 'character',
                    help = 'Protein ID for which the colocalisation is done.')
parser$add_argument('--allele_frequency', type = 'numeric', default = 0.01,
                    help = 'Allele frequency threshold.')
parser$add_argument('--i2', type = 'numeric', default = 40,
                    help = 'I2 threshold.')
parser$add_argument('--sample_fraction_overlap', type = 'numeric', default = 0.5,
                    help = 'Fraction of samples each variants needs to be tested in, out of maximum samples.')

args <- parser$parse_args()

# read in allele data
alleles <- arrow::open_dataset(args$allele_info) %>% select("ID", "CHR", "bp", "str_allele1", "str_allele2")
message("Allele info read in!")
nrow(alleles)

# Read in pQTLs
readpqtl <- function(fn){
    dt_temp <- fread(fn, select = c(3, 1, 2, 4, 5, 10, 11, 8, 6))
    dt_temp <- dt_temp[A1FREQ > args$allele_frequency & A1FREQ < 100 - args$allele_frequency]
    keycols <- c("ID")
    setkeyv(dt_temp, keycols)
    return(dt_temp)
}

pqtl_path <- list.files(args$pqtl_folder, full.names = TRUE)
pqtl_path <- pqtl_path[!str_detect(pqtl_path, "chrX")]
pqtls <- lapply(pqtl_path, readpqtl)
pqtls <- rbindlist(pqtls)

message("pQTLs read in!")

# Read in eQTL data
ds <- arrow::open_dataset(args$eqtl_folder, partitioning = "phenotype", hive_style = TRUE)
eqtls <- ds %>% filter(i_squared < args$i2 & phenotype %in% !!args$gene_id ) %>% collect() %>% 
filter(sample_size >= args$sample_fraction_overlap * max(sample_size)) %>% as.data.table()
message("eQTLs read in!")

alleles <- alleles %>% filter(ID %in% !!eqtls$variant) 
message("Allele info filtered!")
alleles <- alleles %>% collect() %>% as.data.table()
message("Allele info converted!")

eqtls <- merge(eqtls, alleles, by.x = "variant", by.y = "ID")
rm(alleles)
message("Data merged!")

print(head(eqtls))
print(head(pqtls))

# Harmonise datasets
harmonised_data <- harmonise_sumstats(eqtls, 
pqtls, 
data1_chr = "CHR",
data2_chr = "CHROM",
data1_bp = "bp",
data2_bp = "GENPOS",
data1_ea = "str_allele2",
data2_ea = "ALLELE1",
data1_nea = "str_allele1",
data2_nea = "ALLELE0",
data1_beta = "beta",
data2_beta = "BETA",
data1_se = "standard_error",
data2_se = "SE")

saveRDS(harmonised_data, file = paste0(args$protein_id, "__", args$gene_id, ".rds"))
