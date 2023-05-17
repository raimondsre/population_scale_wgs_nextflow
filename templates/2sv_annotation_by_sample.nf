#!/usr/bin/env nextflow
// Pipeline to count structural variants (SVs)
// Counts by sample
// By segment

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
 set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi"), env(variantsPresent) into separated_by_segment

 script:
 input = remExt(vcf.name) 
 """
       bcftools view ${vcf} ${chr}:${start}-${stop} |
       bcftools view --exclude 'POS<${start}' |
       bcftools view --exclude 'POS>${stop}' |
       bcftools view -c1 -Oz -o ${input}.${intervalname}.vcf.gz
       bcftools index -t ${input}.${intervalname}.vcf.gz
       
       variantsPresent=1
       if [ `bcftools view ${input}.${intervalname}.vcf.gz --no-header | wc -l` -eq 0 ]; then variantsPresent=0; fi
 """
}
separated_by_segment = separated_by_segment.filter { it[5] == "1" }.map { tuple(it[0..4]) }
// Separate segment into samples
( separated_by_segment, separated_by_segment_split_samples ) = separated_by_segment.into(2)
separated_by_segment_split_samples = separated_by_segment_split_samples.combine(samples_ch1)
process separateVCF_by_samples {
 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx), val(order_samp), val(sample) from separated_by_segment_split_samples

 output:
 set val(order), val(intervalname), val(input), file("${input}.${intervalname}.${sample}.vcf.gz"), file("${input}.${intervalname}.${sample}.vcf.gz.tbi"), val(sample), env(variantsPresent) into separated_by_segment_and_sample
 script:
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
 conda = '/home/raimondsre/.conda/envs/parallel'
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