#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.sampleLocation = './for.trimming'
params.fastqDir = '.'
params.batchName = 'lv_reference_20220722_502samples'
params.trimGaloreContainer = '/mnt/beegfs2/beegfs_large/raimondsre_add2/genome_analysis/trim_galore_0.6.7.sif'
params.hpc_billing_account = 'bmc_1mgenome'

def remPath(String fileName) {
    return fileName.replaceAll(/.*\//,'')
}

process file_transfer {
    //cpus 2
    clusterOptions "-l nodes=wn61:ppn=1 -A ${params.hpc_billing_account}"

    input:
    tuple val(SAMPLE_ID), val(sample_chunk), val(read1_lftp), val(read2_lftp), val(read1_md5sum), val(read2_md5sum)

    output:
    tuple val(SAMPLE_ID), val(sample_chunk), path(read1), path(read2)

    script:
    read1 = remPath(read1_lftp)
    read2 = remPath(read2_lftp)
    """
    lftp -e "set ssl:verify-certificate no; set net:connection-limit 2; get ${read1_lftp} -o ${read1} & get ${read2_lftp} -o ${read2} & wait all; exit"
    
    md5sum1_r1=${read1_md5sum}
    md5sum2_r1=\$(md5sum ${read1} | awk '{ print \$1 }')
    echo "md5sum read 1 NAS:" \$md5sum1_r1
    echo "md5sum read 1 HPC:" \$md5sum2_r1
       echo -e "---"
    md5sum1_r2=${read2_md5sum}
    md5sum2_r2=\$(md5sum ${read2} | awk '{ print \$1 }')
    echo "md5sum read 2 NAS:" \$md5sum1_r2
    echo "md5sum read 2 HPC:" \$md5sum2_r2
    
    if [ -z "\$md5sum1_r1" ] || [[ "\$md5sum1_r1" == "\$md5sum2_r1" && "\$md5sum1_r2" == "\$md5sum2_r2" ]]; then
    echo "Checksums are equal or missing in NAS."
    else
    echo "Checksums doesn't match."
    rm ${read1}
    rm ${read2}
    exit 1
    fi
    """
}

process adaptor_trimming {
       publishDir params.fastqDir, mode: 'copy', overwrite: true, failOnError: true
    cpus 8
    container params.trimGaloreContainer

    input:
    tuple val(SAMPLE_ID), val(sample_chunk), path(read1), path(read2)

    output:
    tuple path(read1_trimmed), path(read2_trimmed)

    script:
    read1_trimmed = read1.toString().replaceAll("_1.f","_1_val_1.f")
    read2_trimmed = read2.toString().replaceAll("_2.f","_2_val_2.f")
    varCal_tsv = "${params.fastqDir}/${params.batchName}_variant_calling.tsv"
    """
    trim_galore --cores 16 --adapter AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
           --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --quality 20 \
           --paired --no_report_file \
           -o . ${read1} ${read2} 
    
    if [ ! -f ${varCal_tsv} ]; then mkdir -p ${params.fastqDir} && touch ${varCal_tsv}; fi
    echo -e "${SAMPLE_ID}\t0\t0\t${SAMPLE_ID}\t${sample_chunk}\t${params.fastqDir}/${read1_trimmed}\t${params.fastqDir}/${read2_trimmed}" >> ${varCal_tsv}
    """
}

// process save_trimmed {
//     publishDir params.fastqDir, mode: 'move', overwrite: true, failOnError: true

//     input:
//     tuple path(read1_trimmed), path(read2_trimmed)

//     output: 
//     tuple file(read1_trimmed), file(read2_trimmed)

//     script:
//     """
//     read1=\$(readlink -f ${read1_trimmed})
//     read2=\$(readlink -f ${read2_trimmed})
//     mv \$read1 ${params.fastqDir}
//     mv \$read2 ${params.fastqDir}
    
//     """
// }

workflow {

    Channel
        .fromPath(params.sampleLocation)
        .splitCsv(header:false, sep:'\t',strip:true)
        .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5]) }
        .filter { SAMPLE_ID, chunk, read1, read2, read1_md5sum, read2_md5sum ->
            def fastq_path1 = read1
            def fastq_filename = fastq_path1.tokenize('/').last().toString().replaceAll("_1.f","_1_val_1.f")
            def file = new File("${params.fastqDir}/${params.batchName}/${fastq_filename}")
            def file_exists = file.exists()
            if (file_exists) println ">>> WARNING: FILE ${read1.tokenize('/').last()} of ${SAMPLE_ID} sample already trimmed"
            !file_exists
        }
        .set { for_lftp }

    for_lftp
        | file_transfer
        | adaptor_trimming
}
