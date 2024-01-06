#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.sampleLocation = './for.trimming'
params.batchDir = '.'
params.batchName = 'lv_reference_20220722_502samples'
params.trimGaloreContainer = '/mnt/beegfs2/beegfs_large/raimondsre_add2/genome_analysis/trim_galore_0.6.7.sif'
params.hpc_billing_account = 'bmc_1mgenome'

def remPath(String fileName) {
    return fileName.replaceAll(/.*\//,'')
}

process file_transfer {
    cpus 1
    clusterOptions "-l nodes=wn61 -A ${params.hpc_billing_account}"

    input:
    tuple val(SAMPLE_ID), val(sample_chunk), val(read1_lftp), val(read2_lftp)

    output:
    tuple val(SAMPLE_ID), val(sample_chunk), path(read1), path(read2)

    script:
    read1 = remPath(read1_lftp)
    read2 = remPath(read2_lftp)
    """
    lftp -e "set ssl:verify-certificate no; set net:connection-limit 2; get ${read1_lftp} -o ${read1} & get ${read2_lftp} -o ${read2} & wait; exit"
    """
}

process adaptor_trimming {
    cpus 1
    container params.trimGaloreContainer

    input:
    tuple val(SAMPLE_ID), val(sample_chunk), path(read1), path(read2)

    output:
    tuple path(read1_trimmed), path(read2_trimmed)

    script:
    read1_trimmed = read1.toString().replaceAll("_1.f","_1_val_1.f")
    read2_trimmed = read2.toString().replaceAll("_2.f","_2_val_2.f")
    varCal_tsv = "${params.batchDir}/${params.batchName}_variant_calling.tsv"
    batchDir = "${params.batchDir}"
    """
    trim_galore --cores 16 --adapter AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
           --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --quality 20 \
           --paired --no_report_file \
           -o . ${read1} ${read2} 
    
    if [ ! -f ${varCal_tsv} ]; then mkdir -p ${batchDir} && touch ${varCal_tsv}; fi
    echo -e "${SAMPLE_ID}\t0\t0\t${SAMPLE_ID}\t${sample_chunk}\t${batchDir}/${read1_trimmed}\t${batchDir}/${read2_trimmed}" >> ${varCal_tsv}
    """
}

process save_trimmed {
    publishDir params.batchDir, mode: 'move', overwrite: true, failOnError: true

    input:
    tuple path(read1_trimmed), path(read2_trimmed)

    output: 
    tuple path(read1_trimmed), path(read2_trimmed)

    script:
    """
    """
}

workflow {

    Channel
        .fromPath(params.sampleLocation)
        .splitCsv(header:false, sep:'\t',strip:true)
        .map { row -> tuple(row[0], row[1], row[2], row[3]) }
        .filter { SAMPLE_ID, chunk, read1, read2 ->
            def fastq_path1 = read1
            def fastq_filename = fastq_path1.tokenize('/').last().toString().replaceAll("_1.f","_1_val_1.f")
            def file = new File("${params.batchDir}/${params.batchName}/${fastq_filename}")
            def file_exists = file.exists()
            if (file_exists) println ">>> WARNING: FILE ${read1.tokenize('/').last()} of ${SAMPLE_ID} sample already trimmed"
            !file_exists
        }
        .set { for_lftp }

    for_lftp
        | file_transfer
        | adaptor_trimming
        | save_trimmed
}
