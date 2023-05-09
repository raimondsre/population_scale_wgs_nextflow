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
 //.filter({it[1].contains('chr14')})
 //.filter({it[4].contains('chr14_1_5000000')})
 .into { intervals1; intervals2 }
// Samples in first and second input VCF
samples_from_first_and_second_merge_file = vcf_first_extractSamples.combine(vcf_second_extractSamples)
//samples_from_first_and_second_merge_file.subscribe {println it}

process extract_overlapping_vcf_samples {
 input:
 tuple file(vcf1), file(idx1), file(vcf2), file(idx2) from samples_from_first_and_second_merge_file
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
params.outputName = remExt(params.firstVCF)+".merged_with."+remExt(params.secondVCF)


// Make single channel for intervals and vcf file to be imputed
vcfIntervals_first = intervals1.combine(vcf_first)
// Make single channel for intervals and vcf file to be used as imputation panel
vcfIntervals_second = intervals2.combine(vcf_second)

vcfIntervals_first_and_second = vcfIntervals_first.mix(vcfIntervals_second)

//vcfIntervals_first_and_second.subscribe { println it }
 
//###
//### Analysis
//###

process separate_segments {
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
       bcftools view --exclude 'POS>${stop}' -Oz -o ${input}.${intervalname}.vcf.gz
       bcftools index -t ${input}.${intervalname}.vcf.gz
 """
}

separated_by_segment_first_and_second = separated_by_segment_first_and_second
       .groupTuple(by:[0,1])

process merge_segments {
       input:
       set val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment_first_and_second
       output:
       tuple val(order), val(params.outputName), file("merged.${intervalname}.vcf.gz"), file("merged.${intervalname}.vcf.gz.tbi"), env(variantsPresent) into merged_ch
       script:
       first = vcf[0]
       sec = vcf[1]

       """
       bcftools merge ${first} ${sec} -Oz -o merged.${intervalname}.vcf.gz
       bcftools index -t merged.${intervalname}.vcf.gz

       variantsPresent=1
       if [ `bcftools view merged.${intervalname}.vcf.gz --no-header | wc -l` -eq 0 ]; then variantsPresent=0; fi
       """
}

merged_ch
       .filter { it[4] == "1" }
       .map { tuple(it[0..3]) }
       .toSortedList({ a,b -> a[1] <=> b[1] })
       .flatten().buffer ( size: 4 )
       .groupTuple(by:1)
       .into {merged_ch_concat}
/*
// Concatanate segments
process concatanate_segments {
 publishDir params.publishDir, mode: 'copy', overwrite: true

 input:
 set val(order), val(name), file(vcf_all), file(idx_all) from merged_ch_concat 
 output:
 set file(output_full), file(output_full+".tbi")
 script:
 output_full = name+".vcf.gz"
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 bcftools concat -f vcfFiles.txt -Oz -o ${output_full}
 bcftools index -t ${output_full}
 """
}
/*
*/