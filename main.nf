// enabling nextflow DSL v2
nextflow.enable.dsl=2

// date spacing results
timestamoName = (new Date().format("yyyy-MM-dd")) + "-${params.codename}"
params.resultsDir = params.outdir + "/" + timestamoName


motd = """
--------------------------------------------------------------------------
quick-rnaseq-nf ($workflow.manifest.version)
--------------------------------------------------------------------------
Name                        : $params.codename
Session ID                  : $workflow.sessionId
--------------------------------------------------------------------------
RNAseq information
--------------------------------------------------------------------------
Samples file                : $params.experiment.samplesheet
Contrasts                   : $params.experiment.contrasts
Results dir                 : $params.resultsDir
Tx reference                : $params.transcriptome.reference
Tx decoys                   : $params.transcriptome.decoys
Fastp args                  : $params.fastp.args
Salmon index args           : $params.salmon.index.args 
Salmon quant libtype        : $params.salmon.quant.libtype
Salmon quant args           : $params.salmon.quant.args
Counts from abundance       : $params.summarize_to_gene.counts_from_abundance
Organism DB                 : $params.summarize_to_gene.organism_db
QC PCA transform            : $params.qc.pca.transform
QC lfc threshold            : $params.qc.ma.lfc_threshold
QC heamap transform         : $params.qc.sample.transform
DGE threshold               : $params.dge.lfc_threshold
DGE FDR treshold            : $params.dge.fdr
GO organism DB              : $params.gene_ontology.organism_db
GO gene ontology gene id    : $params.gene_ontology.gene_id
GO remove Gencode version   : $params.gene_ontology.remove_gencode_version
GO term FDR threshold       : $params.gene_ontology.fdr
--------------------------------------------------------------------------
Environment information
--------------------------------------------------------------------------
Container       : $workflow.container
Config files    : $workflow.configFiles
Project dir     : $workflow.projectDir
Work dir        : $workflow.workDir
Launch dir      : $workflow.launchDir
Command line    : $workflow.commandLine
Repository      : $workflow.repository
CommitID        : $workflow.commitId
Revision        : $workflow.revision
--------------------------------------------------------------------------
"""

log.info motd


process TRIM_READS {

    tag "${sample}"

    publishDir "${params.resultsDir}/qc/reads", pattern: "*.json", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/qc/reads", pattern: "*.html", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(read1), path(read2)
    output:
        tuple val(sample), path("${read1.simpleName}.trimmed.fastq.gz"), path("${read2.simpleName}.trimmed.fastq.gz"), emit: fastq
        path "${sample}.fastp.json", emit: qc
    
    script:
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

    stub:
    """
    touch ${read1.simpleName}.trimmed.fastq.gz
    touch ${read2.simpleName}.trimmed.fastq.gz
    touch ${sample}.fastp.json  
    """
}

process SALMON_INDEX {
    input: 
        path transcriptome
        path decoys

    output: 
        path 'transcriptome_index'

    script:
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

    stub:
    """
    mkdir transcriptome_index
    """
}

process SALMON_QUANT {

    tag "${sample}"

    publishDir "${params.resultsDir}/quantification/", mode: 'copy', overwrite: true

    input: 
        path transcriptome_index
        tuple val(sample), path(read1), path(read2)

    output:
        path("${sample}")

    script:
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

    stub:
    """
    mkdir ${sample}
    """
}

process SUMMARIZE_TO_GENE {

    publishDir "${params.resultsDir}/dataset/", mode: 'copy', overwrite: true

    input: 
        path sample_sheet
        path salmon_results

    output:
        path 'summarized-experiment.rds'
    
    script:
    """
        quick-rnaseq-summarize-to-gene.R \\
            ${sample_sheet} \\
            -c ${params.summarize_to_gene.counts_from_abundance} \\
            -d ${params.summarize_to_gene.organism_db}
    """

    stub:
    """
    touch summarized-experiment.rds
    """
}

process QC_PCA {

    publishDir "${params.resultsDir}/qc/quantification", mode: 'copy', overwrite: true

    input: 
        path sefile

    output:
        path "pca-plot.pdf"
    
    script:
    """
        quick-rnaseq-pca.R ${sefile} pca-plot.pdf --transform=${params.qc.pca.transform}
    """

    stub:
    """
    touch pca-plot.pdf
    """
}

process QC_COOK {

    publishDir "${params.resultsDir}/qc/dge", mode: 'copy', overwrite: true

    input: 
        path sefile

    output:
        path "cook-plot.pdf"
    
    script:
    """
        quick-rnaseq-cook.R ${sefile} cook-plot.pdf
    """

    stub:
    """
    touch cook-plot.pdf
    """
}

process QC_MAPLOT {

    tag "${contrast1}-vs-${contrast2}"

    publishDir "${params.resultsDir}/qc/dge", mode: 'copy', overwrite: true

    input: 
        path sefile
        tuple val(contrast1), val(contrast2)

    output:
        path "maplot-${contrast1}-vs-${contrast2}.pdf"
    
    script:
    """
        quick-rnaseq-ma.R ${sefile} \\
            maplot-${contrast1}-vs-${contrast2}.pdf \\
            --case ${contrast1} \\
            --control ${contrast2} \\
            -l ${params.qc.ma.lfc_threshold}
    """

    stub:
    """
    touch maplot-${contrast1}-vs-${contrast2}.pdf
    """
}

process QC_SAMPLE {

    publishDir "${params.resultsDir}/qc/quantification", mode: 'copy', overwrite: true

    input: 
        path sefile

    output:
        path "sample-clustering.pdf"
    
    script:
    """
    quick-rnaseq-sample-distance.R ${sefile} sample-clustering.pdf --transform=${params.qc.sample.transform}
    """

    stub:
    """
    touch sample-clustering.pdf
    """
}

process QC_REPORT {
    publishDir "${params.resultsDir}/qc/", mode: 'copy', overwrite: true

    input: 
        path reads_json
        path salmon_quant

    output:
        path "multiqc_report.html"
    
    script:
    """
        multiqc .
    """

    stub:
    """
    touch multiqc_report.html
    """
}

process ANALYSIS_DGE {
    
    tag "${contrast1}-vs-${contrast2}"

    publishDir "${params.resultsDir}/analysis/", mode: 'copy', overwrite: true

    input: 
        path sefile
        tuple val(contrast1), val(contrast2)

    output:
        tuple path("dge-${contrast1}-vs-${contrast2}.csv"), val(contrast1), val(contrast2)
    
    script:
    """
        quick-rnaseq-dge.R \\
            ${sefile} \\
            dge-${contrast1}-vs-${contrast2}.csv \\
            --case ${contrast1} \\
            --control ${contrast2} \\
            -l ${params.dge.lfc_threshold} \\
            -f ${params.dge.fdr}
    """

    stub:
    """
    touch dge-${contrast1}-vs-${contrast2}.csv
    """
}

process ANALYSIS_GO {

    tag "${contrast1}-vs-${contrast2}"

    publishDir "${params.resultsDir}/analysis/", mode: 'copy', overwrite: true

    input: 
        tuple path(results), val(contrast1), val(contrast2)

    output:
        path "go-${contrast1}-vs-${contrast2}.csv"
    
    script:
    """
        quick-rnaseq-go.R ${results} \\
            go-${contrast1}-vs-${contrast2}.csv \\
            -d ${params.gene_ontology.organism_db} \\
            -g ${params.gene_ontology.gene_id} \\
            -f ${params.gene_ontology.fdr} \\
            --remove-gencode-version=${params.gene_ontology.remove_gencode_version}
    """

    stub:
    """
    touch go-${contrast1}-vs-${contrast2}.csv
    """
}

process TELEMETRY {
    publishDir "${params.resultsDir}", mode: 'copy'

    output:
    path('run.info.txt')

    script:
    """
    echo '${motd}' > run.info.txt
    """

    stub:
    """
    touch run.info.txt
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
    QC_COOK(SUMMARIZE_TO_GENE.out)
    QC_REPORT(TRIM_READS.out.qc.collect(), SALMON_QUANT.out.collect())

    // differential expression
    ANALYSIS_DGE(SUMMARIZE_TO_GENE.out, contrasts_ch)

    // gene ontology analysis
    ANALYSIS_GO(ANALYSIS_DGE.out)

    // telemetry information
    TELEMETRY()
}

workflow {
    QUICK_RNASEQ()
}

