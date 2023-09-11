#!/usr/bin/env Rscript

library(SumTool)
library(vroom)
library(rtracklayer)
library(IGUtilityPackage)
library(stringr)
library(argparse)

setDTthreads(1)

parser <- ArgumentParser(description = 'Run LDSC genetic correlation between two datasets.')
parser$add_argument('--harmonised_data', metavar = 'file', type = 'character',
                    help = 'R data file which contains list with harmonised sumstats matrices.')
parser$add_argument('--gene_id', type = 'character',
                    help = 'Gene name for which the LDSC is done.')
parser$add_argument('--protein_id', type = 'character',
                    help = 'Protein ID for which the LDSC is done.')
parser$add_argument('--hapmap3_ld_scores', type = 'character',
                    help = 'Folder with HapMap3 LD-scores.')
parser$add_argument('--chain', type = 'character',
                    help = 'chain file.')

args <- parser$parse_args()

# Function
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
                                    data2_ldscore = "ldcore",
                                    data1_se = "SE",
                                    data1_n = "N",
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
    ldscore = data2[[data2_ldscore]],
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


  print(head(data1))
  print(head(data2))



  if(isTRUE(showProgress)){message("Alleles and effect sizes harmonised!")}
  data1 <- data1[, c(11, 1:7), with = FALSE]
  data2 <- data2[, c(9, 1:5), with = FALSE]
  if(isTRUE(showProgress)){message("Finished!")}
  if(isTRUE(showProgress)){return(list(sumstats1 = data1, sumstats2 = data2))}

}


# Functions
RemoveSigEffects <- function(data, window = 1000000, Pthresh = 5e-8, 
snp_id_col = "SNP", snp_chr_col = "Chr",
snp_pos_col = "Pos", snp_a1_col = "A1", snp_a2_col = "A2",
beta_col = "beta", se_col = "se", n_col = "N") {

  data <- data.table(SNP = data[[snp_id_col]],
                     Chr = data[[snp_chr_col]],
                     Pos = data[[snp_pos_col]],
                     A1 = data[[snp_a1_col]],
                     A2 = data[[snp_a2_col]],
                     BETA = data[[beta_col]],
                     SE = data[[se_col]],
                     N = data[[n_col]]
                     )


  data$P <- ZtoP(data$BETA/data$SE)
  data$Z <- data$BETA/data$SE

  print(head(data))

  #data_f <- data[data$P < as.numeric(Pthresh), ]
  data_f <- data
  # Iteratively identify most significant SNP, and remove all other SNPs in the window
  res <- data_f[-c(1:nrow(data_f)), ]

  while (min(data_f$P) <= 5e-8) {
    lead_snp <- data_f[abs(data_f$Z) == max(abs(data_f$Z)), ]
    if (nrow(lead_snp) > 1) {
      lead_snp <- lead_snp[1, ]
    }
    res <- rbind(res, lead_snp)
    data_f <- data_f[!(data_f$Chr == lead_snp$Chr & data_f$Pos > lead_snp$Pos - window & data_f$Pos < lead_snp$Pos + window), ]
    #message(paste("Added:", lead_snp$snp_chr, lead_snp$snp_pos))
  }
  return(data_f)
}



and <- readRDS(args$harmonised_data)

eqtl <- and$sumstats1
pqtl <- and$sumstats2

colnames(eqtl) <- c("SNP", "Chr", "Pos", "A1", "A2", "BETA", "SE", "N")
colnames(pqtl) <- c("SNP", "Chr", "Pos", "A1", "A2", "BETA", "SE", "N")

# Prepare LD scores
files <- fs::dir_ls(glob = "*ldscore.gz", path = paste0(args$hapmap3_ld_scores, "/", "eur_w_ld_chr"))
ldscore <- vroom(files)
h3 <- vroom(paste0(args$hapmap3_ld_scores, "/", "w_hm3.snplist"))
ldscore <- merge(ldscore, h3)

gr <- GRanges(
  seqnames = Rle(paste0("chr", ldscore$CHR)),
  ranges = IRanges(start = ldscore$BP, end = ldscore$BP),
  strand = Rle(strand("+")),
  SNP = ldscore$SNP,
  A1 = ldscore$A1,
  A2 = ldscore$A2,
  CM = ldscore$CM,
  MAF = ldscore$MAF,
  LD2 = ldscore$L2)

ch = import.chain(args$chain)

gr2 <- liftOver(gr, ch)
gr2 <- as.data.table(unlist(gr2))

ldscore <- data.table(SNP = gr2$SNP,
                      Chr  = str_replace(gr2$seqnames, "chr", ""),
                      Pos = gr2$start,
                      A1 = gr2$A1,
                      A2 = gr2$A2,
                      Maf = gr2$MAF,
                      ldscore = gr2$LD2)

ldscore <- ldscore[order(Chr, Pos)]
ldscore <- ldscore[!duplicated(SNP), ]

message("LD scores lifted")

comb1 <- harmonise_sumstats_fast(data1 = eqtl,
                                data2 = ldscore,
                                data1_chr = "Chr", data2_chr = "Chr",
                                data1_bp = "Pos", data2_bp = "Pos",
                                data1_ea = "A1", data2_ea = "A1",
                                data1_nea = "A2", data2_nea = "A2",
                                data1_beta = "BETA", data2_ldscore = "ldscore",
                                data1_se = "SE", data1_n = "N")

message("1. harmonisation done")

comb2 <- harmonise_sumstats_fast(data1 = pqtl,
                                 data2 = ldscore,
                                 data1_chr = "Chr", data2_chr = "Chr",
                                 data1_bp = "Pos", data2_bp = "Pos",
                                 data1_ea = "A1", data2_ea = "A1",
                                 data1_nea = "A2", data2_nea = "A2",
                                 data1_beta = "BETA", data2_ldscore = "ldscore",
                                 data1_se = "SE", data1_n = "N")

message("2. harmonisation done")

ldscore <- data.table(SNP = comb1$sumstats2$SnpId,
Chr  = comb1$sumstats2$chr,
Pos = comb1$sumstats2$pos,
A1 = comb1$sumstats2$ea,
A2 = comb1$sumstats2$nea,
Maf = 0.1,
ldscore = comb1$sumstats2$ldscore
)

eqtl <- data.table(SNP = comb1$sumstats1$SnpId,
                      Chr  = comb1$sumstats1$chr,
                      Pos = comb1$sumstats1$pos,
                      A1 = comb1$sumstats1$ea,
                      A2 = comb1$sumstats1$nea,
                      BETA = comb1$sumstats1$beta,
                      SE = comb1$sumstats1$se,
                      N = comb1$sumstats1$n)
pqtl <- data.table(SNP = comb2$sumstats1$SnpId,
                   Chr  = comb2$sumstats1$chr,
                   Pos = comb2$sumstats1$pos,
                   A1 = comb2$sumstats1$ea,
                   A2 = comb2$sumstats1$nea,
                   BETA = comb2$sumstats1$beta,
                   SE = comb2$sumstats1$se,
                   N = comb2$sumstats1$n)


gen_cor <- LDreg(sumstat = list(eqtl, pqtl),
              ldscore =  ldscore)

message("Genetic correlation for all loci done.")

h2_1_2 <- gen_cor[3,1]
se_h2_1_2 <- gen_cor[3,2]

h2_1 <- gen_cor[1,1]
h2_2 <- gen_cor[2, 1]
se_h2_1 <- gen_cor[1,1]
se_h2_2 <- gen_cor[2, 1]

rg <- h2_1_2 / (sqrt(h2_1) * sqrt(h2_2))
rg_se <- sqrt((se_h2_1_2 / (sqrt(h2_1) * sqrt(h2_2)))^2 + ((h2_1_2 * se_h2_1 / (2 * h2_1 * sqrt(h2_1) * sqrt(h2_2)))^2 + (h2_1_2 * se_h2_2 / (2 * h2_2 * sqrt(h2_1) * sqrt(h2_2)))^2))

res1 <- data.table(gene = args$gene_id, protein = args$protein_id, trait1_h2 = gen_cor[1, 1], trait1_h2_se = gen_cor[1, 2],
                  trait2_h2 = gen_cor[2, 1], trait2_h2_se = gen_cor[2, 2],
                  trait1_intercept = gen_cor[1, 3], trait1_intercept_se = gen_cor[1, 4],
                  trait2_intercept = gen_cor[2, 3], trait2_intercept_se = gen_cor[2, 4],
                  genetic_covariance = gen_cor[3, 1], genetic_covariance_se = gen_cor[3, 2],
                  trait1_trait2_intercept = gen_cor[3, 3], trait1_trait2_intercept_se = gen_cor[3, 4],
                  rg = rg,
                  rg_se = rg_se,
                  rg_p = ZtoP(rg/rg_se),
                  analysis = "All effects"
                  )
# Recalculate while removing all significant regions

print(head(eqtl))

eqtl2 <- RemoveSigEffects(eqtl, snp_id_col = "SNP", snp_chr_col = "Chr", snp_pos_col = "Pos", beta_col = "BETA", se_col = "SE", n_col = "N", snp_a1_col = "A1", snp_a2_col = "A2")
pqtl2 <- RemoveSigEffects(pqtl, snp_id_col = "SNP", snp_chr_col = "Chr", snp_pos_col = "Pos", beta_col = "BETA", se_col = "SE", n_col = "N", snp_a1_col = "A1", snp_a2_col = "A2")
message("Removal of singificant loci completed.")


eqtl2 <- eqtl2[SNP %in% pqtl$SNP]
pqtl2 <- pqtl2[SNP %in% eqtl$SNP]

gen_cor <- LDreg(sumstat = list(eqtl2[, c(1:8), with = FALSE], pqtl2[, c(1:8), with = FALSE]),
                 ldscore =  ldscore)
message("Genetic correlation for non-sig loci done.")

h2_1_2 <- gen_cor[3,1]
se_h2_1_2 <- gen_cor[3,2]

h2_1 <- gen_cor[1,1]
h2_2 <- gen_cor[2, 1]
se_h2_1 <- gen_cor[1,1]
se_h2_2 <- gen_cor[2, 1]

rg <- h2_1_2 / (sqrt(h2_1) * sqrt(h2_2))
rg_se <- sqrt((se_h2_1_2 / (sqrt(h2_1) * sqrt(h2_2)))^2 + ((h2_1_2 * se_h2_1 / (2 * h2_1 * sqrt(h2_1) * sqrt(h2_2)))^2 + (h2_1_2 * se_h2_2 / (2 * h2_2 * sqrt(h2_1) * sqrt(h2_2)))^2))

res2 <- data.table(gene = args$gene_id, protein = args$protein_id, trait1_h2 = gen_cor[1, 1], trait1_h2_se = gen_cor[1, 2],
                  trait2_h2 = gen_cor[2, 1], trait2_h2_se = gen_cor[2, 2],
                  trait1_intercept = gen_cor[1, 3], trait1_intercept_se = gen_cor[1, 4],
                  trait2_intercept = gen_cor[2, 3], trait2_intercept_se = gen_cor[2, 4],
                  genetic_covariance = gen_cor[3, 1], genetic_covariance_se = gen_cor[3, 2],
                  trait1_trait2_intercept = gen_cor[3, 3], trait1_trait2_intercept_se = gen_cor[3, 4],
                  rg = rg,
                  rg_se = rg_se,
                  rg_p = ZtoP(rg/rg_se),
                  analysis = "Sig. effects removed"
                  )

res <- rbind(res1, res2)

fwrite(res, paste0(args$gene_id, "_", args$protein_id, ".txt"), sep = "\t", quote = FALSE)
message("Analysis finished!")
