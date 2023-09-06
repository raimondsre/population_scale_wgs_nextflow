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
       .groupTuple(by:0)
       .filter {
        a, b, c, d ->
        def number_of_chunks = len(a)
        if (number_of_chunks != 1) println ">>> WARNING: ${a[1]} has multiple chunks, use nf-core/sarek to process"
       }
       .subscribe { println it }
       /*.filter { SAMPLE_ID, chunk, read1, read2 ->
        def number_of_chunks = chunk.toInteger()
        if (number_of_chunks != 1) println ">>> WARNING: ${SAMPLE_ID} has multiple chunks, use nf-core/sarek to process"
        number_of_chunks == 1
       }
       .subscribe { println it }
       //.set { for_lftp }
/*
*/
//test.subscribe { println it }
