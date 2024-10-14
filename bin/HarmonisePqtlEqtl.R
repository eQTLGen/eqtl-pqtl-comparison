#!/usr/bin/env Rscript

# Load libraries
library(IGUtilityPackage)
library(argparse)
library(data.table)
library(stringr)

# functions

harmonise_sumstats_fast <- function(data1, data2,
                                    data1_chr = "CHR",
                                    data2_chr = "CHR",
                                    data1_bp = "BP",
                                    data2_bp = "BP",
                                    data1_ea = "EA",
                                    data2_ea = "NEA",
                                    data1_nea = "EA",
                                    data2_nea = "NEA",
                                    data1_beta = "BETA",
                                    data2_beta = "BETA",
                                    data1_se = "SE",
                                    data2_se = "SE",
                                    data1_n = "N",
                                    data2_n = "N",
                                    showProgress = TRUE){
  data1 <- data.table(
    chr = data1[[data1_chr]],
    pos = data1[[data1_bp]],
    ea = data1[[data1_ea]],
    nea = data1[[data1_nea]],
    beta = data1[[data1_beta]],
    se = data1[[data1_se]],
    n = data1[[data1_n]],
    key = c("chr", "pos")
  )

  data2 <- data.table(
    chr = data2[[data2_chr]],
    pos = data2[[data2_bp]],
    ea = data2[[data2_ea]],
    nea = data2[[data2_nea]],
    beta = data2[[data2_beta]],
    se = data2[[data2_se]],
    n = data2[[data2_n]],
    key = c("chr", "pos")
  )

  if(isTRUE(showProgress)){message("Datasets prepared!")}

  common_pos <- intersect(data1$pos, data2$pos)
  data1 <- data1[pos %in% common_pos]
  data2 <- data2[pos %in% common_pos]

  data1$chr_pos <- paste0(data1$chr, ":", data1$pos)
  data2$chr_pos <- paste0(data2$chr, ":", data2$pos)

  common_pos <- intersect(data1$chr_pos, data2$chr_pos)
  data1 <- data1[chr_pos %in% common_pos]
  data2 <- data2[chr_pos %in% common_pos]

  if(isTRUE(showProgress)){message("Data prefiltered!")}

  # Remove tri-allelic variants
  data1 <- data1[, .SD[!duplicated(chr_pos)]]
  data2 <- data2[, .SD[!duplicated(chr_pos)]]

  if(isTRUE(showProgress)){message("Triallelic variants removed!")}
  # first dataset SNP ID variants
  data1$SNPID1 <- paste0(data1$chr_pos, "_", data1$ea, "_", data1$nea)
  if(isTRUE(showProgress)){message("1 done")}
  data1$SNPID2 <- paste0(data1$chr_pos, "_", data1$nea, "_", data1$ea)
  if(isTRUE(showProgress)){message("2 done")}
  dataset1_variants <- c(data1$SNPID1, data1$SNPID2)
  if(isTRUE(showProgress)){message("combined")}
  # second ID variant
  data2$SNPID1 <- paste0(data2$chr_pos, "_", data2$ea, "_", data2$nea)
  if(isTRUE(showProgress)){message("3 done")}
  data2$SNPID2 <- paste0(data2$chr_pos, "_", data2$nea, "_", data2$ea)
  if(isTRUE(showProgress)){message("4 done")}
  dataset2_variants <- c(data2$SNPID1, data2$SNPID2)
  if(isTRUE(showProgress)){message("combined")}
  # filter
  setkey(data1, SNPID1)
  if(isTRUE(showProgress)){message("Dataset 1 keys set!")}
  setkey(data2, SNPID1)
  if(isTRUE(showProgress)){message("Dataset 2 keys set!")}

  batch_size <- 100000
  num_batches <- ceiling(length(dataset2_variants) / batch_size)
  batches <- split(dataset2_variants, rep(1:num_batches, each = batch_size, length.out = length(dataset2_variants)))

  # Filter the data.table using each batch
  data1 <- rbindlist(lapply(batches, function(batch){data1[SNPID1 %in% batch]}))
  if(isTRUE(showProgress)){message("Dataset 1 filtered!")}

  num_batches <- ceiling(length(dataset1_variants) / batch_size)
  batches <- split(dataset1_variants, rep(1:num_batches, each = batch_size, length.out = length(dataset1_variants)))
  data2 <- rbindlist(lapply(batches, function(batch){data2[SNPID1 %in% batch]}))
  if(isTRUE(showProgress)){message("Dataset 2 filtered!")}

  if(isTRUE(showProgress)){message("Dataset variants overlapped!")}

  # harmonise SNP names
  data1$SnpId <- data1$SNPID1
  data2$SnpId <- data2$SNPID1

  data1[data1$SNPID1 == data2$SNPID2]$beta <- -data1[data1$SNPID1 == data2$SNPID2]$beta
  data1[data1$SNPID1 == data2$SNPID2]$ea <- data2[data1$SNPID1 == data2$SNPID2]$ea
  data1[data1$SNPID1 == data2$SNPID2]$nea <- data2[data1$SNPID1 == data2$SNPID2]$nea
  data1[data1$SNPID1 == data2$SNPID2]$SnpId <- data2[data1$SNPID1 == data2$SNPID2]$SnpId

  if(isTRUE(showProgress)){message("Alleles and effect sizes harmonised!")}
  data1 <- data1[, c(11, 1:7), with = FALSE]
  data2 <- data2[, c(11, 1:7), with = FALSE]
  if(isTRUE(showProgress)){message("Finished!")}
  if(isTRUE(showProgress)){return(list(sumstats1 = data1, sumstats2 = data2))}


}


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
parser$add_argument('--i2', type = 'numeric', default = 100,
                    help = 'I2 threshold.')
parser$add_argument('--sample_fraction_overlap', type = 'numeric', default = 0.5,
                    help = 'Fraction of samples each variants needs to be tested in, out of maximum samples.')

args <- parser$parse_args()

# read in allele data
alleles <- arrow::open_dataset(args$allele_info) %>% select("variant_index", "chromosome", "bp", "non_eff_allele", "eff_allele")
message("Allele info read in!")
nrow(alleles)

# Read in pQTLs
readpqtl <- function(fn){
    dt_temp <- fread(fn, select = c(3, 1, 2, 4, 5, 10, 11, 8, 6))
    dt_temp <- dt_temp[A1FREQ > args$allele_frequency & A1FREQ < 1 - args$allele_frequency]
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
eqtls <- ds %>% filter((i_squared <= args$i2 | is.na(i_squared)) & phenotype %in% !!args$gene_id ) %>% collect() %>% 
filter(sample_size >= args$sample_fraction_overlap * max(sample_size)) %>% as.data.table()
message("eQTLs read in!")

alleles <- alleles %>% filter(variant_index %in% !!eqtls$variant) 
message("Allele info filtered!")
alleles <- alleles %>% collect() %>% as.data.table()
message("Allele info converted!")

eqtls <- merge(eqtls, alleles, by = "variant_index")
rm(alleles)
message("Data merged!")

print(head(eqtls))
print(head(pqtls))

# Harmonise datasets
harmonised_data <- harmonise_sumstats_fast(eqtls, 
pqtls, 
data1_chr = "chromosome",
data2_chr = "CHROM",
data1_bp = "bp",
data2_bp = "GENPOS",
data1_ea = "eff_allele",
data2_ea = "ALLELE1",
data1_nea = "non_ref_allele",
data2_nea = "ALLELE0",
data1_beta = "beta",
data2_beta = "BETA",
data1_se = "standard_error",
data2_se = "SE",
data1_n = "sample_size",
data2_n = "N")

saveRDS(harmonised_data, file = paste0(args$protein_id, "__", args$gene_id, ".rds"))
