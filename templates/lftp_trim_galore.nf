#!/usr/bin/env nextflow

params.sampleLocation = './for.trimming'
params.batchDir = '.'
params.batchName = 'lv_reference_20220722_502samples'

// Read the input file containing sample ID and it's location
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
       //.filter({it[1].contains('25408')})
       .set { for_lftp }
def remPath(String fileName) {return fileName.replaceAll(/.*\//,'')}

for_lftp.subscribe{ println it }


process file_transfer { 
       cpus 1
       clusterOptions 'local'

       input:
       set val(SAMPLE_ID), (sample_chunk), path(read1_lftp), path(read2_lftp) from for_lftp
       
       output:
       set val(SAMPLE_ID), (sample_chunk), file(read1), file(read2) into for_trimming

       script:
       read1 = remPath(${read1_lftp})
       read2 = remPath(${read2_lftp})
       """
       lftp -e "set ssl:verify-certificate no; get ${read1_lftp} -o ${read1}; get ${read2_lftp} -o ${read2}; exit"
       """
}
/*
process adaptor_trimming {
       publishDir params.batchDir+"/"+params.batchName, mode: 'move', overwrite: false
       cpus 16
       
       input:
       set val(SAMPLE_ID), (sample_chunk), path(read1), path(read2) from for_trimming

       output:
       set file(read1_trimmed), file(read2_trimmed)

       script:
       read1_trimmed = read1.toString().replaceAll("_1.f","_1_val_1.f")
       read2_trimmed = read2.toString().replaceAll("_2.f","_2_val_2.f")
       varCal_tsv = "${params.batchDir}/${params.batchName}/${params.batchName}_variant_calling.tsv"
       batchDir = "${params.batchDir}/${params.batchName}"
       """
       trim_galore --cores 16 --adapter AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
              --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --quality 20 \
              --paired --no_report_file \
              -o . ${read1} ${read2}
       
       if [ ! -f ${varCal_tsv} ]; then mkdir -p ${batchDir} && touch ${varCal_tsv}; fi
       echo -e "${SAMPLE_ID}\t0\t0\t${SAMPLE_ID}\t${sample_chunk}\t${batchDir}/${read1_trimmed}\t${batchDir}/${read2_trimmed}" >> ${varCal_tsv}
       """
}

*/