#!/usr/bin/env nextflow

// Pipeline is useful to catch individual sample issues

// Important notes
// 1. Only variants >50 bp, including SVs and MEIs
// 2. Only variants AC >= 1

// Parameters //
params.publishDir = '/Users/work/Documents/ngs/hs/paper202208/pipeline/results'
params.VCFfile = './merged.two.vcf.gz'
params.input = './hg38chr25int5e6.bed'
// Constants
params.baseDir='/home_beegfs/groups/bmc/genome_analysis_tmp/hs'
params.beegfProgs = '/home_beegfs/raimondsre/programmas'
params.annotsvDir = '/home_beegfs/raimondsre/programmas/AnnotSV/bin/AnnotSV'

// Channels //
// Define channels for intervals and initial .vcf.gz file
// VCF input
Channel
 .fromPath(params.VCFfile)
 .map { tuple(it, it+".tbi") }
 .into { vcf; vcf_extractSamples }
// VCF input samples
process extract_vcf_samples {
 input:
 tuple file(vcf), file(idx) from vcf_extractSamples
 output:
 file 'samples' into Channel_samples mode flatten
 script:
 """
 bcftools query -l ${vcf} > samples
 """ }
counter2 = 0
Channel_samples
 .splitText() {it.replaceFirst(/\n/,'')}
 .map {value ->
        counter2 += 1
        [counter2, value].flatten()}
 .into { samples_ch1; samples_ch2}


// Intervals bed
counter = 0
Channel
 .fromPath(params.input)
 .splitCsv(header:false, sep:'\t',strip:true)
 .map { row -> tuple(row[0], row[1], row[2], row[0]+"_"+row[1]+"_"+row[2]) }
 .map {value ->
        counter += 1
        [counter, value].flatten()}
 //.filter({it[1].contains('chrY')})
 .into { intervals1; intervals2 }
// Methods //
// Define function to remove .vcf.gz extension
def remExt(String fileName) {return fileName.replaceFirst(/\.vcf\.gz$/,'')}
// Make single channel for intervals and vcf file
Channel_vcfIntervals = intervals1.combine(vcf)

//###
//### Analysis
//###

// Separate VCF into fragments, has to be before separating out samples
process separateVCF {
 
 input:
 tuple val(order), val(chr), val(start), val(stop), val(intervalname), file(vcf), file(idx) from Channel_vcfIntervals
 
 output:
 set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi") into separated_by_segment

 script:
 input = remExt(vcf.name) 
 """
       bcftools view ${vcf} ${chr}:${start}-${stop} |
       bcftools view --exclude 'POS<${start}' |
       bcftools view --exclude 'POS>${stop}' -Oz -o ${input}.${intervalname}.vcf.gz
       bcftools index -t ${input}.${intervalname}.vcf.gz
 """
}
// Separate segment into samples
( separated_by_segment, separated_by_segment_split_samples ) = separated_by_segment.into(2)
( separated_by_segment_split_samples ) = separated_by_segment_split_samples.combine(samples_ch1)
process separateVCF_by_samples {
 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx), val(order_samp), val(sample) from separated_by_segment_split_samples

 output:
 set val(order), val(intervalname), val(input), file("${input}.${intervalname}.${sample}.vcf.gz"), file("${input}.${intervalname}.${sample}.vcf.gz.tbi"), val(order_samp), val(sample) into separated_by_segment_and_sample
 script:
 """
 bcftools view ${vcf} -s ${sample} -Oz -o ${input}.${intervalname}.${sample}.vcf.gz
 bcftools index -t ${input}.${intervalname}.${sample}.vcf.gz
 """
}

// Annotsv >50 bp variant annotation and counting
process annotsv_sv_annotate_segments {
 
 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment

 output:
 file '*vepCounts.txt' into annotsvCounts_separated_by_segment_ch
 //set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.setID.vcf.gz"), file("${remExt(vcf.name)}.setID.vcf.gz.tbi") into segments_ready_for_collection

 """
 # Variant annotation
 bcftools view -c1 ${vcf} | bcftools sort -Oz -o ${input}.${intervalname}.ac1.vcf.gz
 ${params.annotsvDir} -SVinputFile ${input}.${intervalname}.ac1.vcf.gz \
                -outputFile ${intervalname}.${input}.ac1.annotsv \
                -genomeBuild GRCh38 \
                -overlap 95
  # Count variants by type, function and location
    Rscript annotsv_annotated_variant_counting.R ${intervalname} ${input}
 """
}
process annotsv_sv_annotate_segments_by_sample {
 publishDir = params.publishDir
 
 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx), val(order_samp), val(sample) from separated_by_segment_and_sample

 output:
 file '*vepCounts.txt' into annotsvCounts_separated_by_segment_and_sample_ch
 //set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.setID.vcf.gz"), file("${remExt(vcf.name)}.setID.vcf.gz.tbi"), val(order_samp), val(sample) into segments_sample_ready_for_collection

 """
 chmod ${params.annotsvDir}
  # Variant annotation
 bcftools view -c1 ${vcf} | bcftools sort -Oz -o ${input}.${intervalname}.${sample}.ac1.vcf.gz
 ${params.annotsvDir} -SVinputFile ${input}.${intervalname}.${sample}.ac1.vcf.gz \
                -outputFile ${intervalname}.${input}.ac1.annotsv \
                -genomeBuild GRCh38 \
                -overlap 95
  # Count variants by type, function and location
    Rscript annotsv_annotated_variant_counting.R ${intervalname} ${input}
 """
}

//###
//### Merging
//###

process annotsv_sv_annotate_merge {
       publishDir = params.publishDir
       scratch = '/scratch'
       input:
       file (vep_annotated_counts) from annotsvCounts_separated_by_segment_ch.concat(annotsvCounts_separated_by_segment_and_sample_ch)

       output:
       file 'annotsvCounts.sv.total.txt'
       script:
       """
       cat *annotsvCounts.txt > annotsvCounts.sv.total.txt
       """
}
/*
segments_sample_ready_for_collection_collected = segments_sample_ready_for_collection
 .toSortedList({ a,b -> a[5] <=> b[5] })
 .flatten().buffer( size: 7 )
 .groupTuple(by:[0,1,2]) 

// Merge samples
process merge_samples {
 input:
 set val(order), val(intervalname), val(input), file(vcf_all), file(idx_all), val(order_samp), val(sample_all) from segments_sample_ready_for_collection_collected
 output:
  set val(order), val(intervalname), val(input), file("merged.${intervalname}.vcf.gz"), file("merged.${intervalname}.vcf.gz.tbi"), val(order_samp), val(sample_all) into segments_sample_ready_for_collection_merged

 script:
 """
 bcftools merge ${vcf_all.join(' ')} -Oz -o merged.${intervalname}.vcf.gz
 bcftools index -t merged.${intervalname}.vcf.gz
 """
}

segments_sample_ready_for_collection_collected = segments_sample_ready_for_collection_merged
 .map { tuple(it[0],it[1],it[2],it[3],it[4]) }
 .toSortedList({ a,b -> a[0] <=> b[0] })
 .flatten().buffer ( size: 5 )
 .groupTuple(by:[2])
 
// Arrange segments and group by input file name
segments_ready_for_collection_collected = segments_ready_for_collection
 .toSortedList({ a,b -> a[0] <=> b[0] })
 .flatten().buffer ( size: 5 )
 .groupTuple(by:[2])

// Merge segments
process concatanate_segments {
 publishDir = params.publishDir

 input:
 set val(order), val(intervalname), val(input), file(vcf_all), file(idx_all) from segments_sample_ready_for_collection_collected 
 output:
 file ("merged.vcf.gz")
 script:
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 bcftools concat --naive -f vcfFiles.txt -Oz -o merged.vcf.gz
 """
}
*/
