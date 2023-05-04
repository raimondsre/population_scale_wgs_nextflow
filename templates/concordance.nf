#!/usr/bin/env nextflow
// Testin
params.publishDir = './results'
params.refDir = '/home_beegfs/groups/bmc/genome_analysis_tmp/hs/ref'

params.firstVCF = './'
params.secondVCF = './'
//params.VCFfile = './merged.two.vcf.gz'
//params.intervalsBed = './hg38chr25int5e6.bed'

// Define channels for intervals and initial .vcf.gz file
// Input file for corcondance - first
Channel
 .fromPath(params.firstVCF)
 .map { tuple(it, it+".tbi") }
 .into { vcf_first; vcf_first_extractSamples }
// Input file for corcondance - second
Channel
 .fromPath(params.secondVCF)
 .map { tuple(it, it+".tbi") }
 .into { vcf_second; vcf_second_extractSamples }
// Intervals
counter = 0
Channel
 .fromPath(params.intervalsBed)
 .splitCsv(header:false, sep:'\t',strip:true)
 .map { row -> tuple(row[0], row[1], row[2], row[0]+"_"+row[1]+"_"+row[2]) }
 .map { value ->
        counter += 1
        [counter, value].flatten()}
 .filter({it[1].contains('chrM')})
 .into { intervals1; intervals2 }
// Samples in first and second input VCF
samples_from_first_and_second_concordance_file = vcf_first_extractSamples.combine(vcf_second_extractSamples)
//samples_from_first_and_second_concordance_file.subscribe {println it}

process extract_overlapping_vcf_samples {
 input:
 tuple file(vcf1), file(idx1), file(vcf2), file(idx2) from samples_from_first_and_second_concordance_file
 output:
 file 'samples_overlap' into samples_ch mode flatten
 script:
 """
 bcftools query -l ${vcf1} > samples1
 bcftools query -l ${vcf2} > samples2
 comm -12 <(sort samples1) <(sort samples2) > samples_overlap 
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
def remPath(String fileName) {return fileName.replaceAll(/.*\//,'').replaceFirst(/\.vcf\.gz$/,'')}
def remExt(String fileName) {return fileName.replaceFirst(/\.vcf\.gz$/,'')}
def remExtBref(String fileName) {return fileName.replaceFirst(/\.bref$/,'')}

// Make single channel for intervals and vcf file to be imputed
vcfIntervals_first = intervals1.combine(vcf_first)
// Make single channel for intervals and vcf file to be used as imputation panel
vcfIntervals_second = intervals2.combine(vcf_second)

vcfIntervals_first_and_second = vcfIntervals_first.mix(vcfIntervals_second)

//vcfIntervals_first_and_second.subscribe { println it }
 
//###
//### Analysis
//###

process separateVCF {
 //publishDir params.publishDir

 input:
 tuple val(order), val(chr), val(start), val(stop), val(intervalname), file(vcf), file(idx) from vcfIntervals_first_and_second
 
 output:
 set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi") into separated_by_segment_first_and_second

 script:
 input = remExt(vcf.name) 
 """
       bcftools view ${vcf} ${chr}:${start}-${stop} |
       bcftools view --exclude 'POS<${start}' |
       bcftools view --exclude 'POS>${stop}' |
       bcftools view -v snps -m2 -c3 -Oz -o ${input}.${intervalname}.vcf.gz
       bcftools index -t ${input}.${intervalname}.vcf.gz
 """
}
(separated_by_segment_first_and_second, separated_by_segment_first_and_second_getOverlapID) = separated_by_segment_first_and_second.into(2)

separated_by_segment_first_and_second_getOverlapID = separated_by_segment_first_and_second_getOverlapID
       .groupTuple(by:[0,1])

process finding_overlap_variants {
       input:
       set val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment_first_and_second_getOverlapID
       output:
       tuple val(order), file("variants_overlap.${intervalname}") into overlap_variants
       script:
       first = vcf[0]
       sec = vcf[1]
       """
       bcftools query -f '%ID\n' ${first} > first.id
       bcftools query -f '%ID\n' ${sec} > sec.id
       comm -12 <(sort first.id) <(sort sec.id) > variants_overlap.${intervalname}
       """
}

separated_by_segment_first_and_second_withOverlapID = 
       overlap_variants
       .cross(separated_by_segment_first_and_second)
       .map {tuple(it[1],it[0][1]).flatten()}

// Customise manipulation steps
process manipulate_segment_filtering_overalp_variants {
       //publishDir = params.publishDir
       
       input:
       set val(order), val(intervalname), val(input), file(vcf), file(idx), file(overlap_variants) from separated_by_segment_first_and_second_withOverlapID

       output:
       set val(order), val(intervalname), val(input), file("${input}.${intervalname}.overlapOnly.vcf") into segments_ready_for_concordance

       script:
       """
       bcftools filter -i 'ID=@${overlap_variants}' ${vcf} -Ov -o ${input}.${intervalname}.overlapOnly.vcf
       """
}

segments_ready_for_concordance = segments_ready_for_concordance
       .map { tuple(it, it[2] == remPath(params.firstVCF) ? 0 : 1).flatten() }
       .toSortedList({ a,b -> a[3] <=> b[3] })
       .groupTuple(by:[0,1])
       segments_ready_for_concordance.subscribe {println it}

/*
process manipulate_segment_concordance {
       publishDir = params.publishDir
       
       input:
       set val(order), val(intervalname), val(input), file(vcf) from segments_ready_for_concordance

       output:
       set val(order), val(intervalname), val(input), file("${input}.${intervalname}.overlapOnly.vcf") into segment_concordance

       script:
       """
       SnpSift concordance -v ${vcf[0]} ${vcf[1]} > concordance_${input[0]}_vs_${input[1]}.txt
       """
}
// process manipulate_segment_samples {
//  publishDir = params.publishDir
 
//  input:
//  set val(order), val(intervalname), val(input), file(vcf), file(idx), val(order_samp), val(sample) from separated_by_segment_and_sample

//  output:
//  set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.setID.vcf.gz"), file("${remExt(vcf.name)}.setID.vcf.gz.tbi"), val(order_samp), val(sample) into segments_sample_ready_for_collection

//  """
//  bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${vcf} -Oz -o ${remExt(vcf.name)}.setID.vcf.gz
//  bcftools index -t ${remExt(vcf.name)}.setID.vcf.gz
//  """
//}
/*
//###
//### Merging
//###
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
//segments_sample_ready_for_collection_merged.subscribe {println it}

segments_sample_ready_for_collection_collected = segments_sample_ready_for_collection_merged
 .map { tuple(it[0],it[1],it[2],it[3],it[4]) }
 .toSortedList({ a,b -> a[0] <=> b[0] })
 .flatten().buffer ( size: 5 )
 .groupTuple(by:[2])
//segments_sample_ready_for_collection_collected.subscribe {println it}
 
// Arrange segments and group by input file name
//segments_ready_for_collection_collected = segments_ready_for_collection
// .toSortedList({ a,b -> a[0] <=> b[0] })
// .flatten().buffer ( size: 5 )
// .groupTuple(by:[2])

// Concatanate segments
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