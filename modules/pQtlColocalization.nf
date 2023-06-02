#!/bin/bash nextflow

process HarmoniseSumstats {
    scratch true
    input:
        tuple val(id), val(ensembl), path(pqtl), path(eqtl_folder), path(allele_info_file), path(gtf)

    output:
        tuple val(id), val(ensembl), path(gtf), path("*.rds")

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
    input:
        tuple val(id), val(ensembl), path(gtf), path(harmonised_data)

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
    
    publishDir "${params.outputDir}/ManhattanPlots/", mode: 'copy', overwrite: true

    input:
        tuple val(id), val(ensembl), path(gtf), path(harmonised_data)

    output:
        path "*.png"

    script:
        """
        ManhattanPlot.R \
        --harmonised_data ${harmonised_data} \
        --gene_id ${ensembl} \
        --protein_id ${id} \
        --gtf ${gtf}
        """
}

workflow PQTL_COMPARISON {
    take:
        pqtl_ch

    main:
        harmonised_sumstats_ch = HarmoniseSumstats(pqtl_ch)
        if (params.pQtlEqtlHyrcoloc) {coloc_output_ch = IterativeColoc(harmonised_sumstats_ch)}
        if (params.ManhattanPlots) {manhattan_output_ch = Manhattan(harmonised_sumstats_ch)}

    emit:
        coloc_output_ch
        manhattan_output_ch
}