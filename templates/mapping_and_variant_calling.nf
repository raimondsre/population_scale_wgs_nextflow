#!/usr/bin/env nextflow

// mapping_and_variant_calling.nf
// pipeline is an adaptation of nf-core/sarek pipeline for use on large scale genome analysis
// It performs mapping, aligned read sorting, duplicate marking, raw read quality score recalibration,
// variant calling (haplotypecaller for SNP&INDEL, manta for SV, MELT for MEI)

params.sampleLocation = './for_variant_calling'
Channel
       .fromPath(params.sampleLocation)
       .splitCsv(header:false, sep:'\t',strip:true)
       .map { row -> tuple(row[0], row[1], row[2], row[3]) }
       .subscribe { println }
       /*.filter { SAMPLE_ID, chunk, read1, read2 ->
        if (chunk  1) println ">>> WARNING: ${SAMPLE_ID} has multiple chunks, use nf-core/sarek to process"
        !(chunk > 1)
       }
       .subscribe { println }
       //.set { for_lftp }
for_lftp.subscribe{ println }
/*
*/