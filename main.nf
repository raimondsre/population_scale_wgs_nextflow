#!/usr/bin/env nextflow
// Pipeline to separate VCF by segment, perform custom manipulation and concatanate back
 
params.publishDir = './results'

params.inputVCF = './merged.two.vcf.gz'
params.intervalsBed = './hg38intervals50mil'
params.samplesToKeep = './keep.samples'
params.outputName = remExt(params.inputVCF)+'.filtered'

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
 //.filter({it[1].contains('chrM')})
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
 set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.setID.vcf.gz"), file("${remExt(vcf.name)}.setID.vcf.gz.tbi") into segments_ready_for_collection

 """
 bcftools view -S ${params.samplesToKeep} --force-samples ${vcf} | 
 bcftools view -c3 |
 bcftools norm --multiallelics - |
 bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ${remExt(vcf.name)}.setID.vcf.gz
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
 set file(outputVCF), file(outputVCFtbi) into main_VCF

 script:
 outputVCF = params.outputName+".vcf.gz"
 outputVCFtbi = outputVCF+".tbi"
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 # --naive is risky as it does not check if samples match.
 bcftools concat --naive -f vcfFiles.txt -Oz -o ${outputVCF}
 bcftools index -t ${outputVCF}
 """
}

main_VCF
       .into { main_VCF;
       PCA;
       varinat_counting_all_vep; varinat_counting_all_annotsv;
       varinat_counting_by_sample_vep; varinat_counting_by_sample_annotsv;
       imputation;
       merge_VCF;
       extract_samples
       }

process PCA {
       publishDir params.publishDir, mode: 'move', overwrite: true
       input:
       set file(vcf), file(idx) from PCA
       output:
       file "*.eigenvec" into output_PCA

       shell:
       '''
       nextflow run !{projectDir}/templates/PCA.nf --input !{vcf}
       '''
}

process varinat_counting_all_vep {
       publishDir params.publishDir, mode: 'move', overwrite: true
       input:
       set file(vcf), file(idx) from varinat_counting_all_vep
       output:
       file "*.counted" into output_varinat_counting_all_vep

       shell:
       '''
       nextflow run !{projectDir}/templates/2snp_annotation.nf --input !{vcf}
       '''
}

process varinat_counting_all_annotsv {
       publishDir params.publishDir, mode: 'move', overwrite: true
       input:
       set file(vcf), file(idx) from varinat_counting_all_annotsv
       output:
       file "*.counted" into output_varinat_counting_all_annotsv

       shell:
       '''
       nextflow run !{projectDir}/templates/2sv_annotation.nf --input !{vcf}
       '''
}

process varinat_counting_by_sample_vep {
       publishDir params.publishDir, mode: 'move', overwrite: true
       input:
       set file(vcf), file(idx) from varinat_counting_by_sample_vep
       output:
       file "*.counted" into output_varinat_counting_by_sample_vep

       shell:
       '''
       nextflow run !{projectDir}/templates/2snp_annot_by_sample.nf --input !{vcf}
       '''
}

process varinat_counting_by_sample_annotsv {
       publishDir params.publishDir, mode: 'move', overwrite: true
       input:
       set file(vcf), file(idx) from varinat_counting_by_sample_annotsv
       output:
       file "*.counted" into output_varinat_counting_by_sample_annotsv

       shell:
       '''
       nextflow run !{projectDir}/templates/2sv_annot_by_sample.nf --input !{vcf}
       '''
}

process imputation {
       publishDir params.publishDir, mode: 'move', overwrite: true
       input:
       set file(vcf), file(idx) from imputation
       output:
       file "*.counted.txt" into output_imputation

       shell:
       '''
       nextflow run !{projectDir}/templates/2sv_annot_by_sample.nf --input !{vcf}
       '''
}

process concordance {
       publishDir params.publishDir, mode: 'move', overwrite: true
       input:
       set file(vcf), file(idx) from concordance
       output:
       file "*compared_to*" into output_concordance

       shell:
       '''
       nextflow run !{projectDir}/templates/concordance.nf --input1 !{vcf} --input2 ${params.array_for_concordance}
       '''
}