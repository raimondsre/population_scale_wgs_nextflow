#!/usr/bin/env nextflow
// Pipeline to separate VCF by segment, perform custom manipulation and concatanate back

params.publishDir = './results'

params.inputVCF = './merged.two.vcf.gz'
params.intervalsBed = './hg38intervals50mil'
params.samplesToKeep = './keep.samples'
params.variantsToKeep = './keep.variants'
params.subset = 'all' 
params.outputName = remPath(params.inputVCF)+'.subset_of_'+params.subset
params.phasedDir = '/mnt/beegfs2/home/groups/bmc/references/populationVCF/phased' // Contains phased and bref corrected segments
params.chain_file = '/home/raimondsre/analysis/hs/genome/prs/app/continental_ethnicity/hg38ToHg19.over.chain'
params.target_reference_genome = '/home/raimondsre/analysis/hs/genome/prs/app/continental_ethnicity/Homo_sapiens_assembly38.fasta'

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
 .filter({it[1].contains('22')})
 //.filter({it[4].contains('chr7_155000001_159345973')}) // Imputation problematic with the following segments: chr9_40000001_45000000,
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
def remPath(String fileName) {return fileName.replaceAll(/.*\//,'').replaceFirst(/\.vcf\.gz$/,'')}

// Make single channel for intervals and vcf file
vcfIntervals = intervals1.combine(vcf)
//samples_ch1.subscribe { println it }
 
//###
//### Analysis
//###

// Separate VCF into fragments, has to be before separating by sample
process separateVCF {
 //publishDir params.publishDir
 publishDir params.phasedDir, mode: 'symlink', overwrite: false
 

 input:
 tuple val(order), val(chr), val(start), val(stop), val(intervalname), file(vcf), file(idx) from vcfIntervals
 
 output:
 set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi"), env(variantsPresent) into separated_by_segment

 script:
 input = remExt(vcf.name) 
 """
        if [ -e ${params.phasedDir}/${input}.${intervalname}.vcf.gz ]; then
              ln -s ${params.phasedDir}/${input}.${intervalname}.vcf.gz ${input}.${intervalname}.vcf.gz
              ln -s ${params.phasedDir}/${input}.${intervalname}.vcf.gz.tbi ${input}.${intervalname}.vcf.gz.tbi
       else
              bcftools view ${vcf} ${chr}:${start}-${stop} |
              bcftools view --exclude 'POS<${start}' |
              bcftools view --exclude 'POS>${stop}' -Oz -o ${input}.${intervalname}.vcf.gz
              bcftools index -t ${input}.${intervalname}.vcf.gz
       fi



       if [ "${remPath(params.chain_file)}" = "hg19ToHg38.over.chain.gz" ] || [ "${remPath(params.chain_file)}" = "hg19ToHg38.over.chain" ]; then
       for CHR in {1..22} X Y MT; do
          echo ${CHR} chr${CHR}
          done >> chr_names.txt
       bcftools annotate --rename-chrs chr_names.txt ${input}.${intervalname}.vcf.gz -Oz -o chromosome_corrected.${input}.${intervalname}.vcf.gz
       mv chromosome_corrected.${input}.${intervalname}.vcf.gz ${input}.${intervalname}.vcf.gz
       bcftools index -t ${input}.${intervalname}.vcf.gz
       fi

       variantsPresent=1
       if [ `bcftools view ${input}.${intervalname}.vcf.gz --no-header | wc -l` -eq 0 ]; then variantsPresent=0; fi
 """
}
separated_by_segment 
       .filter { it[5] == "1" }
       .into {separated_by_segment_filtered}
// Customise manipulation steps
process manipulate_segment {
 //publishDir params.publishDir
 cpus 4
 time 24.h
 memory 16.GB
 container '/home/raimondsre/analysis/hs/genome/prs/app/continental_ethnicity/picard:3.1.1--hdfd78af_0'
 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment_filtered

 output:
 set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.setID.vcf.gz"), file("${remExt(vcf.name)}.setID.vcf.gz.tbi") into segments_ready_for_collection

 """
       export _JAVA_OPTIONS="-Xmx16g"
       picard LiftoverVcf \
       I=${vcf} \
       O=${remExt(vcf.name)}.setID.vcf.gz \
       WARN_ON_MISSING_CONTIG=true \
       CHAIN=${params.chain_file} \
       REJECT=${remExt(vcf.name)}.setID.rejected_variants.vcf.gz \
       R=${params.target_reference_genome}
 """
}
// bcftools +setGT ${vcf} -- -t q -n . -i 'FMT/GQ<20' |
// bcftools +setGT ${vcf} -- -t q -n . -i 'FMT/DP<30' |
// bcftools view -i 'F_MISSING <= 0.1' | 
// bcftools +fill-tags -Oz -o ${remExt(vcf.name)}.setID.vcf.gz -- -t ExcHet,AC,AF
 
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
 outputVCF = params.outputName+".vcf.gz"
 outputVCFtbi = outputVCF+".tbi"
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 # --naive is risky as it does not check if samples match.
 bcftools concat -a -f vcfFiles.txt -Oz -o ${outputVCF}
 bcftools index -t ${outputVCF}
 """
}

