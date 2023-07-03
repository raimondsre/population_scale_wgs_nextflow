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
       .into { for_trimming }
def remPath(String fileName) {return fileName.replaceAll(/.*\//,'')}


// Encrypt each file
process adaptor_trimming {
       publishDir params.batchDir+"/"+params.batchName, mode: 'move', overwrite: false
       conda 'trim'
       cpus 16
       
       input:
       set val(SAMPLE_ID), (sample_chunk), path(read1), path(read2) from for_trimming

       output:
       set file(read1_trimmed), file(read2_trimmed) into read_encrypted

       script:
       read1_trimmed = read1.toString().replaceAll("_1.","_1_val_1.")
       read2_trimmed = read2.toString().replaceAll("_2.","_2_val_2.")
       """
       trim_galore --cores 16 --adapter AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
                --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --quality 20 \
                --paired --no_report_file \
                -o . ${read1} ${read2}
       """
}