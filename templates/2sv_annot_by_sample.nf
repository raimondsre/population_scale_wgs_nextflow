#!/usr/bin/env nextflow
// main_bySegmentOnly.nf
params.publishDir = './results'
params.annotsvDir = '/home_beegfs/raimondsre/programmas/AnnotSV'

params.inputVCF = './merged.two.vcf.gz'
params.intervalsBed = './intervals50mil'
// Define channels for intervals and initial .vcf.gz file
// Input file
Channel
 .fromPath(params.inputVCF)
 .map { tuple(it, it+".tbi") }
 .into { vcf; vcf_extractSamples }
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

//###
//### Analysis
//###

// Separate segment into samples
separated_by_segment_split_samples = vcf
 .map { tuple(1, "all_intervals", it[0], it[1]) }
 .combine(samples_ch1)

separated_by_segment_split_samples.subscribe {println it}
/*
process separateVCF_by_samples {
 input:
 set val(order), val(intervalname), file(vcf), file(idx), val(order_samp), val(sample) from separated_by_segment_split_samples

 output:
 set val(order), val(intervalname), val(input), file("${input}.${intervalname}.${sample}.vcf.gz"), file("${input}.${intervalname}.${sample}.vcf.gz.tbi"), val(sample), env(variantsPresent) into separated_by_segment_and_sample
 script:
 input = remExt(vcf.name)
 """
 bcftools view ${vcf} -s ${sample} |
 bcftools view -c1 -Oz -o ${input}.${intervalname}.${sample}.vcf.gz
 bcftools index -t ${input}.${intervalname}.${sample}.vcf.gz

 variantsPresent=1
 if [ `bcftools view ${input}.${intervalname}.${sample}.vcf.gz --no-header | wc -l` -eq 0 ]; then variantsPresent=0; fi
 """
}
separated_by_segment_and_sample = separated_by_segment_and_sample.filter { it[6] == "1" }.map { tuple(it[0..5]) }


// add segment channel sameple element and concatanate with segments_samples channel
//separated_by_segment_and_sample = separated_by_segment.map { tuple(it[0..4],"all").flatten() }.mix(separated_by_segment_and_sample)

// Customise manipulation steps
process manipulate_segment_by_interval_and_sample_annotsv {
 //publishDir params.publishDir
 //cpus 1
 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx), val(sample) from separated_by_segment_and_sample

 output:
 set val(order), val(intervalname), val(input), file("${intervalname}.${sample}.annotsv.counted") into segments_ready_for_collection

 script:
 vcf_name = vcf.name
 """
 uname -a | awk '{print \$2}'

 export ANNOTSV=${params.annotsvDir}
 ${params.annotsvDir}/bin/AnnotSV -SVinputFile ${vcf} \
                -outputFile ${intervalname}.${input}.ac1.annotsv \
                -outputDir . \
                -genomeBuild GRCh38 \
                -overlap 95
 # Count features
 Rscript ${projectDir}/countANNOTSVfeatures.R --input ${intervalname}.${input}.ac1.annotsv.tsv --interval ${intervalname} --sample ${sample} --original_file_name ${input}
 """
}

// Arrange segments and group by input file name
segments_ready_for_collection_collected = segments_ready_for_collection
 //sorts by the first variable [0]
 .toSortedList({ a,b -> a[0] <=> b[0] })
 //flattens whole channel and then remakes sets by 4 elements (if input sets were 5 elements, number should be changed)
 .flatten().buffer ( size: 4 )
 //squishes all elements into one set, collects. For each set variable [2] there will be one channel element.
 .groupTuple(by:[1,2])

// Concatanate segments
process concatanate_segments_by_interval {
 //publishDir params.publishDir, mode: 'move', overwrite: true
 //cpus 16
 input:
 set val(order), val(intervalname), val(input), file(annotsv_all) from segments_ready_for_collection_collected 
 output:
 set val(input), file("${input}.${intervalname}.annotsv.counted") into concatanate_segments_whole
 script:
 """
 cat ${annotsv_all.join(' ')} > ${input}.${intervalname}.annotsv.counted
 """
}

concatanate_segments_whole = concatanate_segments_whole.groupTuple(by:[0])
process concatanate_segments {
 publishDir params.publishDir, mode: 'move', overwrite: true
 //cpus 16
 input:
 set val(input), file(annotsv_all) from concatanate_segments_whole 
 output:
 set file("${input}.by_sample.annotsv.counted")
 script:
 """
 cat ${annotsv_all.join(' ')} > ${input}.by_sample.annotsv.counted
 """
}