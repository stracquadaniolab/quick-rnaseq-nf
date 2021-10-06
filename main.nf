// enabling nextflow DSL v2
nextflow.enable.dsl=2

process TRIM_READS {

    tag "${sample}"

    publishDir "${params.resultsdir}/qc/reads", pattern: "*.json", mode: 'copy', overwrite: true
    publishDir "${params.resultsdir}/qc/reads", pattern: "*.html", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(read1), path(read2)
    output:
        tuple val(sample), path("${read1.simpleName}.trimmed.fastq.gz"), path("${read2.simpleName}.trimmed.fastq.gz"), emit: fastq
        path "${sample}.fastp.json", emit: qc
    
    """
        fastp -w ${task.cpus} \\
            ${params.fastp.args}  \\
            --in1 ${read1} \\
            --in2 ${read2} \\
            --out1 ${read1.simpleName}.trimmed.fastq.gz \\
            --out2 ${read2.simpleName}.trimmed.fastq.gz \\
            --json ${sample}.fastp.json \\
            --html ${sample}.fastp.html
    """

}

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

    publishDir "${params.resultsdir}/quantification/", mode: 'copy', overwrite: true

    input: 
        path transcriptome_index
        tuple val(sample), path(read1), path(read2)

    output:
        path("${sample}")

    """
        salmon quant \\
        -p ${task.cpus} \\
        -l ${params.salmon.quant.libtype} \\
        -i ${transcriptome_index} \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o ${sample} \\
        ${params.salmon.quant.args}
    """
}

process SUMMARIZE_TO_GENE {

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

    publishDir "${params.resultsdir}/qc/quantification", mode: 'copy', overwrite: true

    input: 
        path sefile

    output:
        path "pcaplot.pdf"
    
    """
        quick-rnaseq-pca.R ${sefile} pcaplot.pdf --transform=${params.qc.pca.transform}
    """
}

process QC_MAPLOT {

    tag "${contrast1}-vs-${contrast2}"

    publishDir "${params.resultsdir}/qc/dge", mode: 'copy', overwrite: true

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

process QC_SAMPLE {

    publishDir "${params.resultsdir}/qc/quantification", mode: 'copy', overwrite: true

    input: 
        path sefile

    output:
        path "sample-clustering.pdf"
    
    """
        quick-rnaseq-sample-distance.R ${sefile} sample-clustering.pdf --transform=${params.qc.sample.transform}
    """
}

process QC_REPORT {
    publishDir "${params.resultsdir}/qc/", mode: 'copy', overwrite: true

    input: 
        path reads_json
        path salmon_quant

    output:
        path "multiqc_report.html"
    
    """
        multiqc .
    """
}

process ANALYSIS_DGE {
    
    tag "${contrast1}-vs-${contrast2}"

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

    tag "${contrast1}-vs-${contrast2}"

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
    
    // sample channels
    samplesheet_file = file(params.experiment.samplesheet)
    samples_ch = channel.from(samplesheet_file)
                    .splitCsv(header: true)
                    .map{ record -> tuple(record.sample, file(record.read1), file(record.read2)) }

    // reads trimming
    TRIM_READS(samples_ch)

    // quantification
    SALMON_QUANT(SALMON_INDEX.out, TRIM_READS.out.fastq)

    // gene level quantification
    SUMMARIZE_TO_GENE(samplesheet_file, SALMON_QUANT.out.collect())

    // deseq qc
    contrasts_ch = channel.from(params.experiment.contrasts).map{ x,y -> tuple(x,y) }
    QC_PCA(SUMMARIZE_TO_GENE.out)
    QC_SAMPLE(SUMMARIZE_TO_GENE.out)
    QC_MAPLOT(SUMMARIZE_TO_GENE.out, contrasts_ch)
    QC_REPORT(TRIM_READS.out.qc.collect(), SALMON_QUANT.out.collect())

    // differential expression
    ANALYSIS_DGE(SUMMARIZE_TO_GENE.out, contrasts_ch)

    // gene ontology analysis
    ANALYSIS_GO(ANALYSIS_DGE.out)
}

workflow {
    QUICK_RNASEQ()
}

