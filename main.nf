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

    publishDir "${params.resultsDir}/salmon/", pattern: "${sample}_quant", mode: 'copy', overwrite: true

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

    publishDir "${params.resultsDir}/deseq/", mode: 'copy', overwrite: true

    input: 
        path sample_sheet
        path salmon_results

    output:
        path 'summarized-experiment.rds'
    
    """
        summarize_to_gene.R 
    """
}

workflow {
    // transcriptome indexing
    tx = file(params.transcriptome)
    SALMON_INDEX(tx)
    
    // mRna quantification
    sample_sheet_file = file(params.sampleSheet)
    samples_ch = Channel.from(sample_sheet_file)
                    .splitCsv(header: true)
                    .map{ record -> tuple(record.sample, file(record.read1), file(record.read2)) }

    SALMON_QUANT(SALMON_INDEX.out, samples)

    // gene level quantification
    SUMMARIZE_TO_GENE(sample_sheet_file, SALMON_QUANT.out.collect())
}

