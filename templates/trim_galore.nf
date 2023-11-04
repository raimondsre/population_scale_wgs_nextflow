#!/usr/bin/env nextflow
// prepares input for nf-core/sarek pipeline

params.sampleLocation = './for.trimming'
params.batchDir = '.'
//params.batchName = 'lv_reference_20220722_502samples'

// Read the input file containing sample ID and it's location
Channel
       .fromPath(params.sampleLocation)
       .splitCsv(header:false, sep:'\t',strip:true)
       .map { row -> tuple(row[0], row[1]) }
       //.filter({it[1].contains('25408')})
       .set { for_trimming }
def remPath(String fileName) {return fileName.replaceAll(/.*\//,'')}

//for_trimming.subscribe{ println it }



process adaptor_trimming {
       publishDir params.batchDir, mode: 'move', overwrite: false
       cpus 16
       executor = 'pbs'
       
       input:
       set path(read1), path(read2) from for_trimming

       output:
       set file(read1_trimmed), file(read2_trimmed)

       script:
       read1_trimmed = read1.toString().replaceAll("_1.f","_1_val_1.f")
       read2_trimmed = read2.toString().replaceAll("_2.f","_2_val_2.f")
       """
       trim_galore --cores 16 --adapter AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
              --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --quality 20 \
              --paired --no_report_file \
              -o . ${read1} ${read2}
       """
}
