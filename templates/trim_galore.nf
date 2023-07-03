#!/usr/bin/env nextflow

params.sampleLocation = './for.trimming'
params.batchDir = '.'
params.batchName = 'lv_reference_20220722_502samples'

// Read the input file containing sample ID and it's location
Channel
       .fromPath(params.sampleLocation)
       .splitCsv(header:false, sep:'\t',strip:true)
       .map { row -> tuple(row[0], row[1], row[2], row[3]) }
       //.filter({it[1].contains('25408')})
       .set { for_trimming }
def remPath(String fileName) {return fileName.replaceAll(/.*\//,'')}


// Encrypt each file
process adaptor_trimming {
       publishDir params.batchDir+"/"+params.batchName, mode: 'move', overwrite: false
       cpus 16
       
       input:
       set val(SAMPLE_ID), (sample_chunk), path(read1), path(read2) from for_trimming

       output:
       set file(read1_trimmed), file(read2_trimmed) into read_encrypted

       script:
       read1_trimmed = read1.toString().replaceAll("_1.f","_1_val_1.f")
       read2_trimmed = read2.toString().replaceAll("_2.f","_2_val_2.f")
       """
       trim_galore --cores 16 --adapter AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
                --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --quality 20 \
                --paired --no_report_file \
                -o . ${read1} ${read2}
       echo -e "${SAMPLE_ID}\t0\t0${SAMPLE_ID}\t${sample_chunk}\t${params.batchDir}/${params.batchName}/${read1_trimmed}\t${params.batchDir}/${params.batchName}/${read1_trimmed}" >> ${params.batchDir}/${params.batchName}/variant_calling.tsv
       """
}