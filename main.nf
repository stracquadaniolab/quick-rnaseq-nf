// enabling nextflow DSL v2
nextflow.enable.dsl=2

process SALMON_INDEX {
    input: 
        path transcriptome
        path decoys

    output: 
        path 'transcriptome_index'

    """ 
        # build decoys id list
        grep '^>' <(gunzip -c ${decoys}) | cut -d " " -f 1 > decoys.txt
        sed -i.bak -e 's/>//g' decoys.txt

        # build gentrome
        cat ${transcriptome} ${decoys} > gentrome.fa.gz

        salmon index \\
            -p ${task.cpus} \\
            -t gentrome.fa.gz \\
            -d decoys.txt \\
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
        -l ${params.salmon.quant.libtype} \\
        -i ${transcriptome_index} \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o ${sample}_quant \\
        ${params.salmon.quant.args}
    """
}

process SUMMARIZE_TO_GENE {

    label 'process_high'

    publishDir "${params.resultsdir}/dataset/", mode: 'copy', overwrite: true

    input: 
        path sample_sheet
        path salmon_results

    output:
        path 'summarized-experiment.rds'
    
    """
        quick-rnaseq-summarize-to-gene.R \\
            ${sample_sheet} \\
            -c ${params.summarize_to_gene.counts_from_abundance} \\
            -d ${params.summarize_to_gene.organism_db}
    """
}

process QC_PCA {

    publishDir "${params.resultsdir}/qc/", mode: 'copy', overwrite: true

    input: 
        path sefile

    output:
        path "pcaplot.pdf"
    
    """
        quick-rnaseq-pca.R ${sefile} pcaplot.pdf
    """
}

process QC_MAPLOT {

    publishDir "${params.resultsdir}/qc/", mode: 'copy', overwrite: true

    input: 
        path sefile
        tuple val(contrast1), val(contrast2)

    output:
        path "maplot-${contrast1}-${contrast2}.pdf"
    
    """
        quick-rnaseq-ma.R ${sefile} \\
            maplot-${contrast1}-${contrast2}.pdf \\
            --case ${contrast1} \\
            --control ${contrast2} 
    """
}


process ANALYSIS_DGE {

    publishDir "${params.resultsdir}/analysis/", mode: 'copy', overwrite: true

    input: 
        path sefile
        tuple val(contrast1), val(contrast2)

    output:
        tuple path("dexp-${contrast1}-${contrast2}.csv"), val(contrast1), val(contrast2)
    
    """
        quick-rnaseq-dge.R \\
            ${sefile} \\
            dexp-${contrast1}-${contrast2}.csv \\
            --case ${contrast1} \\
            --control ${contrast2}             
    """
}

process ANALYSIS_GO {

    publishDir "${params.resultsdir}/analysis/", mode: 'copy', overwrite: true

    input: 
        tuple path(results), val(contrast1), val(contrast2)

    output:
        path "go-${contrast1}-${contrast2}.csv"
    
    """
        quick-rnaseq-go.R \\
            ${results} \\ 
            go-${contrast1}-${contrast2}.csv \\
            -d ${params.gene_ontology.organism_db} \\
            -g ${params.gene_ontology.gene_id} \\ 
    """
}

workflow QUICK_RNASEQ{
    
    // decoy-aware transcriptome indexing
    tx_ref_file = file(params.transcriptome.reference)
    tx_decoy_file = file(params.transcriptome.decoys)
    SALMON_INDEX(tx_ref_file, tx_decoy_file)
    
    // mRNA quantification
    samplesheet_file = file(params.experiment.samplesheet)
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
    ANALYSIS_DGE(SUMMARIZE_TO_GENE.out, contrasts_ch)

    // gene ontology analysis
    ANALYSIS_GO(ANALYSIS_DGE.out)
}

workflow {
    QUICK_RNASEQ()
}

