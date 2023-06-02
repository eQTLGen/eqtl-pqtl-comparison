#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2


def helpmessage() {

log.info"""

HyprColocForEqtlCatalog v${workflow.manifest.version}"
===========================================================
Pipeline for running HyprColoc colocalisation analyses (https://www.nature.com/articles/s41467-020-20885-8) between GWAS locus and eQTL Catalogue datasets (https://www.ebi.ac.uk/eqtl/).

Usage:

nextflow run RunHyprColocOnGWAS.nf --eqtl_files \
--pqtl_files \
--allele_info \
--genes \
--pqtl_meta_table \
--gtf  \
--gtf \
--OutputDir

Mandatory arguments:
--eqtl_files              eQTLGen parquet dataset.
--pqtl_files              folder with pQTL files from Sun et. al.
--allele_info             File with allele info for eQTLGen dataset.
--genes                   File with ENSG names to include.
--pqtl_meta_table         Table with protein-gene annotations.
--allele_info             Parquet file with alleles and SNP positions for eQTL dataset.
--gtf                     GTF file for gene annotation.

Optional arguments:
--outputDir               Output directory. Defaults to "results".
--pQtlEqtlHyrcoloc        Run iterative colocalisation. Defaults to true.
--ManhattanPlots          Make Manhattan plots. Defaults to false.
""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.outputDir = 'results'

params.inclusion_step_output = 'NO_FILE'

params.pQtlEqtlHyrcoloc = true
params.ManhattanPlots = false

//Show parameter values
log.info """=======================================================
Colocalisation pipeline v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Output directory']                         = params.OutputDir
summary['eQTL results']                             = params.eqtl_files
summary['pQTL files']                               = params.pqtl_files
summary['Genes']                                    = params.genes
summary['Meta table']                               = params.pqtl_meta_table
summary['Allele info file']                         = params.allele_info
summary['GTF file']                                 = params.gtf
summary['Run hypcoloc']                             = params.pQtlEqtlHyrcoloc
summary['Make Manhattan plots']                     = params.ManhattanPlots

// import modules
include { PQTL_COMPARISON; IterativeColoc; HarmoniseSumstats } from './modules/pQtlColocalization.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
// Get eQTL channel
empirical_ch = Channel.fromPath(params.eqtl_files, type: 'dir')
allelel_ch = Channel.fromPath(params.allele_info)
empirical_ch = empirical_ch.combine(allelel_ch)
genes_ch = Channel.fromPath(params.genes).splitCsv(header: ['gene']).map { row -> tuple(0, row.gene) }
gtf_ch = Channel.fromPath(params.gtf)

// pQTL arguments
pqtl_files_ch = Channel.fromPath(params.pqtl_files).
    map { file ->
          def fileSplit = file.name.toString().split('_')
          def assay = fileSplit[0]
          def uniprot = fileSplit[1]
          def olinkId = fileSplit[2]

          def key = "${assay}_${uniprot}_${olinkId}"

          return tuple(key, file) }

pqtl_ch = Channel.fromPath(params.pqtl_meta_table).splitCsv(header: true, sep: '\t')
    .map { row ->
        def key = "${row.Assay}_${row.UniProt}_${row.OlinkID}"
        return tuple(key, row.ensembl_id)
        }
    .join(genes_ch, by: 1)
    .map { row -> tuple(row[1], row[0])}
    .view()
    .join(pqtl_files_ch)


//empirical_ch.view()
pqtl_ch = pqtl_ch.combine(empirical_ch)
pqtl_ch = pqtl_ch.combine(gtf_ch)

// Define parameter channels
window = Channel.value(params.window)
p_thresh = Channel.value(params.p_thresh)
maf_thresh = Channel.value(params.maf_thresh)
posterior_threshold = Channel.value(params.posterior_threshold)
cs_threshold = Channel.value(params.cs_threshold)

workflow {

        PQTL_COMPARISON(pqtl_ch)
        if (params.pQtlEqtlHyrcoloc) {
        PQTL_COMPARISON.out.coloc_output_ch.flatten().collectFile(name: 'PqtlEqtlColocResults.txt', keepHeader: true, sort: true, storeDir: "${params.outputDir}")
        }

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
