#!/usr/bin/env nextflow
// main_bySegmentOnly.nf
params.publishDir = './results'

params.inputVCF = './merged.two.vcf.gz'
params.intervalsBed = './hg38intervals50mil'
params.samplesToKeep = 'all'
// Define channels for intervals and initial .vcf.gz file
// Input file
Channel
 .fromPath(params.inputVCF)
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
       bcftools view --exclude 'POS<${start}' |
       bcftools view --exclude 'POS>${stop}' -Oz -o ${input}.${intervalname}.vcf.gz
       bcftools index -t ${input}.${intervalname}.vcf.gz
 """
}

// Customise manipulation steps
process manipulate_segment {
 //publishDir params.publishDir
 cpus 1

 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment

 output:
 set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.filtered.vcf.gz"), file("${remExt(vcf.name)}.filtered.vcf.gz.tbi") into segments_ready_for_collection

 """
 if [ ${params.samplesToKeep} != 'all' ]; then 
 bcftools view -S ${params.samplesToKeep} --force-samples ${vcf} -Oz -o ${remExt(vcf.name)}.selected_samples.vcf.gz
 mv ${remExt(vcf.name)}.selected_samples.vcf.gz ${vcf}
 bcftools index -tf ${vcf}
 fi

 bcftools +fill-tags -- -t AF,AC,F_MISSING ${vcf} |
 bcftools view -i 'F_MISSING<0.1' |
 bcftools view -i 'AF>0.01' -Oz -o ${remExt(vcf.name)}.filtered.vcf.gz
 bcftools index -t ${remExt(vcf.name)}.filtered.vcf.gz
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
 set file(outputVCF), file(outputVCFtbi)

 script:
 outputVCF = input+".vcf.gz"
 outputVCFtbi = outputVCF+".tbi"
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 # --naive should to be used with caution.
 bcftools concat --naive -f vcfFiles.txt -Oz -o ${outputVCF}
 bcftools index -t ${outputVCF}
 """
}

