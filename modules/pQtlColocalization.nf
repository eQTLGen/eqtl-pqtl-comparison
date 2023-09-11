#!/bin/bash nextflow

process HarmoniseSumstats {

    container = 'quay.io/urmovosa/pqtlvseqtl:v0.2'
    //scratch true
    input:
        tuple val(id), val(ensembl), path(pqtl), path(eqtl_folder), path(allele_info_file), path(gtf), path(hapmap3_ld_scores), path(chain)

    output:
        tuple val(id), val(ensembl), path(gtf), path("*.rds"), path(hapmap3_ld_scores), path(chain)

    shell:
        '''
        mkdir untar
        tar -C untar -xvf !{pqtl}
        pqtl_path=!{pqtl}
        pqtl_path=$(echo "${pqtl_path%.*}")

        mkdir tmp_eqtls

        cp -r "!{eqtl_folder}/phenotype=!{ensembl}" tmp_eqtls/

        HarmonisePqtlEqtl.R \
        --pqtl_folder untar/${pqtl_path} \
        --eqtl_folder !{eqtl_folder} \
        --allele_info !{allele_info_file} \
        --gene_id !{ensembl} \
        --protein_id !{id}

        rm -r untar
        rm -r tmp_eqtls
        '''
}

process IterativeColoc {

    container = 'quay.io/urmovosa/pqtlvseqtl:v0.2'
    input:
        tuple val(id), val(ensembl), path(gtf), path(harmonised_data), path(hapmap3_ld_scores), path(chain)

    output:
        path "*.txt"

    script:
        """
        IterativeHyprcoloc.R \
        --harmonised_data ${harmonised_data} \
        --gene_id ${ensembl}\
        --protein_id ${id}
        """
}

process Manhattan {
    
    container = 'quay.io/urmovosa/pqtlvseqtl:v0.2'
    publishDir "${params.outdir}/ManhattanPlots/", mode: 'copy', overwrite: true

    input:
        tuple val(id), val(ensembl), path(gtf), path(harmonised_data), path(hapmap3_ld_scores), path(chain)

    output:
        path "*.png", emit: manhattan_output_ch

    script:
        """
        ManhattanPlot.R \
        --harmonised_data ${harmonised_data} \
        --gene_id ${ensembl} \
        --protein_id ${id} \
        --gtf ${gtf}
        """
}


process GeneticCorrelation {

    container = 'quay.io/urmovosa/pqtlvseqtl:v0.2'
    
    input:
        tuple val(id), val(ensembl), path(gtf), path(harmonised_data), path(hapmap3_ld_scores), path(chain)

    output:
        path "*.txt", emit: ldsc_output_ch

    script:
        """
        GeneticCorrelation.R \
        --harmonised_data ${harmonised_data} \
        --gene_id ${ensembl} \
        --protein_id ${id} \
        --hapmap3_ld_scores ${hapmap3_ld_scores} \
        --chain ${chain}
        """
}

process GenomeWideMR {
    
    input:
        tuple val(id), val(ensembl), path(gtf), path(harmonised_data), path(hapmap3_ld_scores), path(chain)

    output:
        path "*.txt", emit: gwmr_output_ch

    script:
        """
        GenomeWideTwoSampleMr.R \
        --harmonised_data ${harmonised_data} \
        --gene_id ${ensembl} \
        --protein_id ${id}
        """
}


workflow HARMONISE {
    take:
        data

    main:
        harmonised_sumstats_ch = HarmoniseSumstats(data)
        
    emit:
        harmonised_sumstats_ch

}

workflow COLOC {
    take:
        data

    main:
        coloc_output_ch = IterativeColoc(data)

    emit:
        coloc_output_ch
}

workflow MANHATTAN {
    take:
        data

    main:
       manhattan_output_ch = Manhattan(data)

    emit:
        manhattan_output_ch
}

workflow LDSC {
    take:
        data

    main:
       ldsc_output_ch = GeneticCorrelation(data)

    emit:
        ldsc_output_ch
}

workflow GWMR {
    take:
        data

    main:
       gwmr_output_ch = GenomeWideMR(data)

    emit:
        gwmr_output_ch
}