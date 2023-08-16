#!/usr/bin/env Rscript

library(hyprcoloc)
library(argparse)
library(data.table)
library(IGUtilityPackage)
#library(Cairo)

parser <- ArgumentParser(description = 'Run HyprColoc for every locus to detect colocalisation between pQTL and eQTL datasets.')
parser$add_argument('--harmonised_data', metavar = 'file', type = 'character',
                    help = 'R data file which contains list with harmonised sumstats matrices')
parser$add_argument('--gene_id', type = 'character',
                    help = 'Gene name for which the colocalisation is done.')
parser$add_argument('--protein_id', type = 'character',
                    help = 'Protein ID for which the colocalisation is done.')
parser$add_argument('--gtf ', type = 'character',
                    help = 'GTF file.')

args <- parser$parse_args()

harmonised_data <- readRDS(args$harmonised_data)

eqtl <- harmonised_data$sumstats1
pqtl <- harmonised_data$sumstats2

eqtl$P <- ZtoP(eqtl$beta/eqtl$se)
pqtl$P <- ZtoP(pqtl$beta/pqtl$se)


# png(paste0(args$protein_id, "__", args$gene_id, "_Manhattan.png"), 
# height = 7, width = 8, res = 400, units = "in", 
# type = "cairo")

# par(mfrow = c(2, 1), mar = c(4, 4, 4, 4))

# ManhPlot(eqtl,
#         InputChrCol = "chr",
#         InputPosCol = "pos",
#         InputPvalCol = "P",
#         InputRsCol = "SnpId",
#         build = "hg38",
#         title = paste("eQTL:", args$gene_id), 
#         AnnotateClosestGenes = TRUE,
#         colors = c("steelblue1", "steelblue3"),
#         gtf = args$gtf)

# ManhPlot(pqtl,
#         InputChrCol = "chr",
#         InputPosCol = "pos",
#         InputPvalCol = "P",
#         InputRsCol = "SnpId",
#         build = "hg38",
#         title = paste("pQTL:", args$protein_id),
#         AnnotateClosestGenes = TRUE,
#         colors = c("palegreen1", "palegreen3"),
#         gtf = args$gtf)

# dev.off()


pdf(paste0(args$protein_id, "__", args$gene_id, "_Manhattan.pdf"), 
height = 7, width = 8)

par(mfrow = c(2, 1), mar = c(4, 4, 4, 4))

ManhPlot(eqtl,
        InputChrCol = "chr",
        InputPosCol = "pos",
        InputPvalCol = "P",
        InputRsCol = "SnpId",
        build = "hg38",
        title = paste("eQTL:", args$gene_id), 
        AnnotateClosestGenes = TRUE,
        colors = c("#00A5FE", "#0072B2"),
        gtf = args$gtf)

ManhPlot(pqtl,
        InputChrCol = "chr",
        InputPosCol = "pos",
        InputPvalCol = "P",
        InputRsCol = "SnpId",
        build = "hg38",
        title = paste("pQTL:", args$protein_id),
        AnnotateClosestGenes = TRUE,
        colors = c("#FB811F", "#D55E00"),
        gtf = args$gtf)

dev.off()