#!/usr/bin/env nextflow
// Pipeline to count variants in VCF file in selected subsets (currently: by INFO (R2) score, by AF bracket)

params.publishDir = './results'
params.refDir = '/home_beegfs/groups/bmc/genome_analysis_tmp/hs/ref'
params.phasedDir = '/mnt/beegfs2/home/groups/bmc/references/populationVCF/phased' // Contains phased and bref corrected segments
params.cpus = 8

params.imputedGenotypes = './'
//params.VCFfile = './merged.two.vcf.gz'
//params.intervalsBed = './hg38chr25int5e6.bed'

// Define channels for intervals and initial .vcf.gz file
// Input file to be imputed
Channel
 .fromPath(params.imputedGenotypes)
 .map { tuple(it, it+".tbi") }
 .into { vcf_imputedGenotypes }
// Intervals
counter = 0
Channel
 .fromPath(params.intervalsBed)
 .splitCsv(header:false, sep:'\t',strip:true)
 .map { row -> tuple(row[0], row[1], row[2], row[0]+"_"+row[1]+"_"+row[2]) }
 .map { value ->
        counter += 1
        [counter, value].flatten()}
 .filter { !(it[1] in ['chrX','chrY','chrM']) }
 //.filter({it[1].contains('chr22')})
 //.filter({it[4].contains('chr4_190000001_190214555')}) // Imputation problematic with the following segments: chr9_40000001_45000000,
 .into { intervals1; intervals2 }
 
// Define function to remove .vcf.gz extension
def remPath(String fileName) {return fileName.replaceAll(/.*\//,'').replaceFirst(/\.vcf\.gz$/,'')}
def remExt(String fileName) {return fileName.replaceFirst(/\.vcf\.gz$/,'')}
def remExtBref(String fileName) {return fileName.replaceFirst(/\.bref$/,'')}

// Make single channel for intervals and vcf file to be imputed
vcfIntervals_imputedGenotypes = intervals1.combine(vcf_imputedGenotypes)

//###
//### Analysis
//###

// Separate VCF into fragments (has to be before separating by sample)
process separateVCF {
       //publishDir params.phasedDir, mode: 'copy', overwrite: false
       input:
       tuple val(order), val(chr), val(start), val(stop), val(intervalname), file(vcf), file(idx) from vcfIntervals_imputedGenotypes
       
       output:
       set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi"), env(variantsPresent) into separated_by_segment

       script:
       input = remExt(vcf.name) 
       """
              bcftools view ${vcf} ${chr}:${start}-${stop} |
              bcftools view --exclude 'POS<${start}' |
              bcftools view --exclude 'POS>${stop}' |
              bcftools norm --remove-duplicates |
              bcftools +fill-tags -- -t AF,AC |
              bcftools annotate --rename-annots <(echo "INFO/R2 INFO") |
              bcftools view -c1 -Oz -o ${input}.${intervalname}.vcf.gz
              bcftools index -t ${input}.${intervalname}.vcf.gz
       variantsPresent=1 
       # Check wether VCF segment has variants present
       if [ `bcftools view ${input}.${intervalname}.vcf.gz --no-header | wc -l` -lt 50 ]; then variantsPresent=0; fi
       """
}
separated_by_segment = separated_by_segment
       .filter { it[5] == "1" }
// Counting variant number by info score
process count_by_info_score {
       //publishDir params.publishDir, mode: 'copy', overwrite: true
       input:
       set val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment

       output:
       set val(name), file(output_counted) into counted_segments_ready_for_collection

       script:
       output_full = remExt(vcf.name)+".txt"
       output_counted = remExt(vcf.name)+".counted.txt"
       name = "${remExt(vcf.name)}" - ".${intervalname}"
       """
       echo -e 'CHR\tSNP\tREF\tALT\tAF\tINFO\tAC' > ${output_full}
       
       bcftools query -f '%CHROM\t%CHROM\\_%POS\\_%REF\\_%ALT\t%REF\t%ALT\t%INFO/AF\t%INFO/INFO\t%INFO/AC\n' ${vcf} >> ${output_full}
       Rscript ${projectDir}/count_info_imputation.R --input ${output_full} --interval ${intervalname} --original_file_name ${remExt(vcf.name)}
       """
}
counted_segments_ready_for_collection = counted_segments_ready_for_collection
       .groupTuple()

process count_by_info_score_collected {
       publishDir params.publishDir, mode: 'copy', overwrite: true
       input:
       set val(name), file(counted_all) from counted_segments_ready_for_collection

       output:
       file output_full

       script:
       output_full = name+".counted.txt"
       """       

       echo -e 'AF_GROUP\tsnv\tINFO_GROUP\tcount\tinterval\tsource' > ${output_full}
       cat ${counted_all.join(' ')} >> ${output_full}
       """
}

/*
*/