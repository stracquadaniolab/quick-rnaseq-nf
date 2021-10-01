// enabling nextflow DSL v2
nextflow.enable.dsl=2

process SALMON_INDEX {
    input: 
        path transcriptome

    output: 
        path 'transcriptome_index'

    """ 
        salmon index \\
            -p ${task.cpus} \\
            -t ${transcriptome} \\
            -i transcriptome_index \\
            ${params.salmon.index.args} 
    """

}

process SALMON_QUANT {

    tag "${sample}"

    label "process_high"

    publishDir "${params.resultsdir}/salmon/", pattern: "${sample}_quant", mode: 'copy', overwrite: true

    input: 
        path transcriptome_index
        tuple val(sample), path(read1), path(read2)

    output:
        path("${sample}_quant")

    """
        salmon quant \\
            -p ${task.cpus} \\
            -l A \\
            -i transcriptome_index \\
            -1 ${read1} \\
            -2 ${read2} \\
            -o ${sample}_quant \\
            --validateMappings
    """
}

process SUMMARIZE_TO_GENE {

    label 'process_high'

    publishDir "${params.resultsdir}/quant/", mode: 'copy', overwrite: true

    input: 
        path sample_sheet
        path salmon_results

    output:
        path 'summarized-experiment.rds'
    
    """
        summarize_to_gene.R ${sample_sheet}
    """
}

process QC_PCA {

    publishDir "${params.resultsdir}/qc/", mode: 'copy', overwrite: true

    input: 
        path sefile

    output:
        path "pcaplot.pdf"
    
    """
        qc_pca.R ${sefile} pcaplot.pdf
    """
}

process QC_MAPLOT {

    publishDir "${params.resultsdir}/qc/", mode: 'copy', overwrite: true

    input: 
        path sefile
        tuple val(ccase), val(ccontrol)

    output:
        path "maplot-${ccase}-${ccontrol}.pdf"
    
    """
        qc_maplot.R ${sefile} maplot-${ccase}-${ccontrol}.pdf --control ${ccontrol} --case ${ccase}
    """
}


process DIFFERENTIAL_EXPRESSION {

    publishDir "${params.resultsdir}/analysis/", mode: 'copy', overwrite: true

    input: 
        path sefile
        tuple val(ccase), val(ccontrol)

    output:
        tuple path("dexp-${ccase}-${ccontrol}.csv"), val(ccase), val(ccontrol)
    
    """
        differential_expression.R ${sefile} dexp-${ccase}-${ccontrol}.csv --control ${ccontrol} --case ${ccase}
    """
}

process GO_ANALYSIS {

    publishDir "${params.resultsdir}/analysis/", mode: 'copy', overwrite: true

    input: 
        tuple path(results), val(ccase), val(ccontrol)

    output:
        path "go-${ccase}-${ccontrol}.csv"
    
    """
        qrna_go.R ${results} go-${ccase}-${ccontrol}.csv
    """
}

workflow {

    // channel.from(params.differential_expression.contrasts).map{ x,y -> tuple(x,y) }.view { x,y -> "$x vs $y"}
    
    // decoy-aware transcriptome indexing
    tx_ref_file = file(params.transcriptome.reference)
    // tx_decoy_file = file(params.transcriptome.decoys)
    SALMON_INDEX(tx_ref_file)
    
    // mRNA quantification
    samplesheet_file = file(params.experiment.samplesheet)
    print(samplesheet_file)
    samples_ch = channel.from(samplesheet_file)
                    .splitCsv(header: true)
                    .map{ record -> tuple(record.sample, file(record.read1), file(record.read2)) }

    SALMON_QUANT(SALMON_INDEX.out, samples_ch)

    // gene level quantification
    SUMMARIZE_TO_GENE(samplesheet_file, SALMON_QUANT.out.collect())

    // deseq qc
    contrasts_ch = channel.from(params.experiment.contrasts).map{ x,y -> tuple(x,y) }
    QC_PCA(SUMMARIZE_TO_GENE.out)
    QC_MAPLOT(SUMMARIZE_TO_GENE.out, contrasts_ch)

    // differential expression
    DIFFERENTIAL_EXPRESSION(SUMMARIZE_TO_GENE.out, contrasts_ch)

    // gene ontology analysis
    GO_ANALYSIS(DIFFERENTIAL_EXPRESSION.out)
}

