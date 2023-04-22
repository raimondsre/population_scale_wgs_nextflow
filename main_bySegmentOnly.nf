#!/usr/bin/env nextflow
// Testin
params.publishDir = './results'

params.VCFfile = './merged.two.vcf.gz'
params.intervalsBed = './hg38intervals50mil'

// Define channels for intervals and initial .vcf.gz file
// Input file
Channel
 .fromPath(params.VCFfile)
 .map { tuple(it, it+".tbi") }
 .into { vcf; vcf_extractSamples }
// Intervals
counter = 0
Channel
 .fromPath(params.intervalsBed)
 .splitCsv(header:false, sep:'\t',strip:true)
 .map { row -> tuple(row[0], row[1], row[2], row[0]+"_"+row[1]+"_"+row[2]) }
 .map {value ->
        counter += 1
        [counter, value].flatten()}
 //.filter({it[1].contains('chrY')})
 .into { intervals1; intervals2 }
// Samples in VCF
process extract_vcf_samples {
 input:
 tuple file(vcf), file(idx) from vcf_extractSamples
 output:
 file 'samples' into samples_ch mode flatten
 script:
 """
 bcftools query -l ${vcf} > samples
 """
}
counter2 = 0
samples_ch
 .splitText() {it.replaceFirst(/\n/,'')}
 .map {value ->
        counter2 += 1
        [counter2, value].flatten()}
 .into { samples_ch1; samples_ch2}

// Define function to remove .vcf.gz extension
def remExt(String fileName) {return fileName.replaceFirst(/\.vcf\.gz$/,'')}

// Make single channel for intervals and vcf file
vcfIntervals = intervals1.combine(vcf)
//samples_ch1.subscribe { println it }
 
//###
//### Analysis
//###

// Separate VCF into fragments, has to be before separating by sample
process separateVCF {
 
 input:
 tuple val(order), val(chr), val(start), val(stop), val(intervalname), file(vcf), file(idx) from vcfIntervals
 
 output:
 set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi") into separated_by_segment

 script:
 input = remExt(vcf.name) 
 """
 bcftools view ${vcf} ${chr}:${start}-${stop} -Oz -o ${input}.${intervalname}.vcf.gz
 bcftools index -t ${input}.${intervalname}.vcf.gz
 """
}

// Separate segment into samples
// ( separated_by_segment, separated_by_segment_split_samples ) = separated_by_segment.into(2)
// separated_by_segment_split_samples = separated_by_segment_split_samples.combine(samples_ch1)
// process separateVCF_by_samples {
//  input:
//  set val(order), val(intervalname), val(input), file(vcf), file(idx), val(order_samp), val(sample) from separated_by_segment_split_samples

//  output:
//  set val(order), val(intervalname), val(input), file("${input}.${intervalname}.${sample}.vcf.gz"), file("${input}.${intervalname}.${sample}.vcf.gz.tbi"), val(order_samp), val(sample) into separated_by_segment_and_sample
//  script:
//  """
//  bcftools view ${vcf} -s ${sample} -Oz -o ${input}.${intervalname}.${sample}.vcf.gz
//  bcftools index -t ${input}.${intervalname}.${sample}.vcf.gz
//  """
// }

// Customise manipulation steps
process manipulate_segment {
 //publishDir = params.publishDir
 
 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment

 output:
 set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.setID.vcf.gz"), file("${remExt(vcf.name)}.setID.vcf.gz.tbi") into segments_ready_for_collection

 """
 bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${vcf} -Oz -o ${remExt(vcf.name)}.setID.vcf.gz
 bcftools index -t ${remExt(vcf.name)}.setID.vcf.gz
 """
}

// process manipulate_segment_samples {
//  //publishDir = params.publishDir
 
//  input:
//  set val(order), val(intervalname), val(input), file(vcf), file(idx), val(order_samp), val(sample) from separated_by_segment_and_sample

//  output:
//  set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.setID.vcf.gz"), file("${remExt(vcf.name)}.setID.vcf.gz.tbi"), val(order_samp), val(sample) into segments_sample_ready_for_collection

//  """
//  bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${vcf} -Oz -o ${remExt(vcf.name)}.setID.vcf.gz
//  bcftools index -t ${remExt(vcf.name)}.setID.vcf.gz
//  """
// }

//###
//### Merging
//###
// segments_sample_ready_for_collection_collected = segments_sample_ready_for_collection
//  .toSortedList({ a,b -> a[5] <=> b[5] })
//  .flatten().buffer( size: 7 )
//  .groupTuple(by:[0,1,2]) 
// // Merge samples
// process merge_samples {
//  input:
//  set val(order), val(intervalname), val(input), file(vcf_all), file(idx_all), val(order_samp), val(sample_all) from segments_sample_ready_for_collection_collected
//  output:
//   set val(order), val(intervalname), val(input), file("merged.${intervalname}.vcf.gz"), file("merged.${intervalname}.vcf.gz.tbi"), val(order_samp), val(sample_all) into segments_sample_ready_for_collection_merged

//  script:
//  """
//  bcftools merge ${vcf_all.join(' ')} -Oz -o merged.${intervalname}.vcf.gz
//  bcftools index -t merged.${intervalname}.vcf.gz
//  """
// }
//segments_sample_ready_for_collection_merged.subscribe {println it}

// segments_sample_ready_for_collection_collected = segments_sample_ready_for_collection_merged
//  .map { tuple(it[0],it[1],it[2],it[3],it[4]) }
//  .toSortedList({ a,b -> a[0] <=> b[0] })
//  .flatten().buffer ( size: 5 )
//  .groupTuple(by:[2])
//segments_sample_ready_for_collection_collected.subscribe {println it}
 
// Arrange segments and group by input file name
segments_ready_for_collection_collected = segments_ready_for_collection
.toSortedList({ a,b -> a[0] <=> b[0] })
.flatten().buffer ( size: 5 )
.groupTuple(by:[2])

// Concatanate segments
process concatanate_segments {
 publishDir params.publishDir, mode: 'move', overwrite: false

 input:
 set val(order), val(intervalname), val(input), file(vcf_all), file(idx_all) from segments_ready_for_collection_collected 
 output:
 file ("merged.vcf.gz")
 script:
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 bcftools concat --naive -f vcfFiles.txt -Oz -o merged.vcf.gz
 """
}
