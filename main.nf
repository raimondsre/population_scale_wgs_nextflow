#!/usr/bin/env nextflow
// main_bySegmentOnly.nf
params.publishDir = './results'

params.VCFfile = './merged.two.vcf.gz'
params.intervalsBed = './hg38intervals50mil'
params.samplesToKeep = './keep.samples'
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
 //publishDir params.publishDir

 input:
 tuple val(order), val(chr), val(start), val(stop), val(intervalname), file(vcf), file(idx) from vcfIntervals
 
 output:
 set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi") into separated_by_segment

 script:
 input = remExt(vcf.name) 
 """
 bcftools view ${vcf} ${chr}:${start}-${stop} |
 bcftools view --exclude 'POS<${start} && POS>${stop}' -Oz -o ${input}.${intervalname}.vcf.gz
 bcftools index -t ${input}.${intervalname}.vcf.gz
 """
}

// Customise manipulation steps
process manipulate_segment {
 //publishDir params.publishDir
 cpus 2

 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment

 output:
 set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.setID.vcf.gz"), file("${remExt(vcf.name)}.setID.vcf.gz.tbi") into segments_ready_for_collection

 """
 bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${vcf} |
 bcftools view -S ${params.samplesToKeep} --force-samples | 
 bcftools norm --multiallelics - |
 bcftools view -c3 -Oz -o ${remExt(vcf.name)}.setID.vcf.gz
 bcftools index -t ${remExt(vcf.name)}.setID.vcf.gz
 """
}

// Arrange segments and group by input file name
segments_ready_for_collection_collected = segments_ready_for_collection
.toSortedList({ a,b -> a[0] <=> b[0] })
.flatten().buffer ( size: 5 )
.groupTuple(by:[2])

// Concatanate segments
process concatanate_segments {
 publishDir params.publishDir, mode: 'move', overwrite: true
 //cpus 16
 input:
 set val(order), val(intervalname), val(input), file(vcf_all), file(idx_all) from segments_ready_for_collection_collected 
 output:
 file ("merged.vcf.gz")
 script:
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 # --naive is risky as it does not check if samples match.
 bcftools concat --naive -f vcfFiles.txt -Oz -o merged.vcf.gz
 """
}
