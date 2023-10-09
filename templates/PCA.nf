#!/usr/bin/env nextflow
// Pipeline to perform PCA analysis of VCF file
// Option --samplesToKeep allows for PCA of subsets

params.publishDir = './results'

params.inputVCF = './merged.two.vcf.gz'
params.intervalsBed = './hg38intervals50mil'
params.samplesToKeep = 'all'
params.variant_missingness_rate = 0.1
params.proportion_of_variants_present = 0.1
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
 .filter { !(it[1] in ['chrX','chrY','chrM']) }
 //.filter({it[1].contains('chr22')})
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
def remPath(String fileName) {return fileName.replaceAll(/.*\//,'').replaceFirst(/\.vcf\.gz$/,'')}
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
       bcftools view --exclude 'POS>${stop}' -Oz -o ${input}.${intervalname}.vcf.gz
       bcftools index -t ${input}.${intervalname}.vcf.gz

       variantsPresent=1
       # Check wether VCF segment has variants present
       if [ `bcftools view ${input}.${intervalname}.vcf.gz --no-header | wc -l` -eq 0 ]; then variantsPresent=0; fi
 """
}
separated_by_segment = separated_by_segment.filter { it[5] == "1"  }.map { tuple(it[0..4]) }

// Customise manipulation steps
process manipulate_segment {
 //publishDir params.publishDir
 cpus 1

 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment

 output:
 set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.filtered.vcf.gz"), file("${remExt(vcf.name)}.filtered.vcf.gz.tbi") into segments_ready_for_collection

 script:
 input = input+".subset_of_"+remPath(params.samplesToKeep)
 """
 if [ ${params.samplesToKeep} != 'all' ]; then 
 bcftools view -S ${params.samplesToKeep} --force-samples ${vcf} -Oz -o ${remExt(vcf.name)}.selected_samples.vcf.gz
 mv ${remExt(vcf.name)}.selected_samples.vcf.gz ${vcf}
 bcftools index -tf ${vcf}
 fi

 bcftools +fill-tags ${vcf} -- -t AF,AC,F_MISSING |
 bcftools view -i 'F_MISSING<${params.variant_missingness_rate}' |
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
 publishDir params.publishDir, mode: 'copy', overwrite: true

 input:
 set val(order), val(intervalname), val(input), file(vcf_all), file(idx_all) from segments_ready_for_collection_collected 
 output:
 set val(input), file(outputVCF) into plink_conversion

 script:
 outputVCF = input+".filtered.vcf.gz"
 outputVCFtbi = outputVCF+".tbi"
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 # --naive should to be used with caution.
 bcftools concat --naive -f vcfFiles.txt -Oz -o ${outputVCF}
 # no indexing necessary for plink2
 # bcftools index -t ${outputVCF}
 """
}

// Convertion to PLINK, additional filtering
process plink_conversion {
       conda = '/home/raimondsre/.conda/envs/plink'
       cpus 32

       input:
       set val(input), file(vcf) from plink_conversion

       output:
       set val(input), file("${input}.{bim,bed,fam}") into pca_analysis

       script:
       """
       plink2 --vcf ${vcf} --threads 32 \
       --output-chr chrM --rm-dup force-first \
       --snps-only --vcf-half-call h --max-alleles 2 \
       --make-bed --out ${input}

       """
}

// PCA analysis
process pca_analysis {
       conda = '/home/raimondsre/.conda/envs/plink'
       publishDir params.publishDir, mode: 'move', overwrite: true
       cpus 8

       input:
       set val(input), file(plink) from pca_analysis

       output:
       file("${input}.eigenvec")

       script:
       """
       plink2 --bfile ${input} --mind ${params.proportion_of_variants_present} --hwe 1e-6 --geno 0.02 --make-bed --out ${input}
       
       plink2 --bfile ${input} --allow-extra-chr --indep-pairwise 50 10 0.5 --out ${input}; 
       plink2 --bfile ${input} --allow-extra-chr --allow-no-sex --chr {1..22} --extract  ${input}.prune.in --pca 10 --out ${input}
      
       """
}
