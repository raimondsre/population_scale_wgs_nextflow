#!/usr/bin/env nextflow
// Testin
params.publishDir = './results'
params.refDir = '/home_beegfs/groups/bmc/genome_analysis_tmp/hs/ref'
params.phasedDir = '/mnt/beegfs2/home/groups/bmc/references/populationVCF/phased'
params.cpus = 32

params.toBeImputed = './'
params.imputationPanel1 = './'
//params.VCFfile = './merged.two.vcf.gz'
//params.intervalsBed = './hg38chr25int5e6.bed'

// Define channels for intervals and initial .vcf.gz file
// Input file to be imputed
Channel
 .fromPath(params.toBeImputed)
 .map { tuple(it, it+".tbi") }
 .into { vcf_toBeImputed; vcf_toBeImputed_extractSamples }
// Input file to be used as imputation panel
Channel
 .fromPath(params.imputationPanel1)
 .map { tuple(it, it+".tbi") }
 .into { vcf_imputation_panel }
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
 .filter({it[4].contains('chr22_50000001_50818468')}) //select the shortest interval of ch18
 .into { intervals1; intervals2 }
 
// Samples in input VCF
process extract_vcf_samples {
 input:
 tuple file(vcf), file(idx) from vcf_toBeImputed_extractSamples
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
def remExtBref(String fileName) {return fileName.replaceFirst(/\.bref$/,'')}

// Make single channel for intervals and vcf file to be imputed
vcfIntervals_toBeImputed = intervals1.combine(vcf_toBeImputed)
// Make single channel for intervals and vcf file to be used as imputation panel
vcfIntervals_toBeUsedAsImputationPanel = intervals2.combine(vcf_imputation_panel)
// Concatanate to-be-separated channels into a single channel
vcfIntervals_toBeImputed_and_toBeUsedAsImputationPanel = vcfIntervals_toBeImputed.mix(vcfIntervals_toBeUsedAsImputationPanel)

//###
//### Analysis
//###

// Separate VCF into fragments (has to be before separating by sample)
process separateVCF {
 //publishDir = params.publishDir
 input:
 tuple val(order), val(chr), val(start), val(stop), val(intervalname), file(vcf), file(idx) from vcfIntervals_toBeImputed_and_toBeUsedAsImputationPanel
 
 output:
 set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi"), env(variantsPresent) into separated_by_segment_toBeImputed_and_toBeUsedAsImputationPanel

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
       # Check wether genetic map file contains variants in segment
       if [ `awk '\$2 >= ${start} && \$2 <= ${stop} {print \$0}' ${params.refDir}/imputation/mapChr/eagle_${chr}_b38.map | wc -l` -eq 0 ]; then variantsPresent=0; fi
 """
}

// Select toBeImputed - impPanel segment pairs only, if both have at least 1 variant after filtering
separated_by_segment_toBeImputed_and_toBeUsedAsImputationPanel
       .filter { it[5] == "1" }
       .groupTuple(by:0)
       .filter { it[2].size() == 2 } //only if both passed
       .multiMap { 
              ch_one: tuple (it[0], it[1][0], it[2][0], it[3][0], it[4][0])
              ch_two: tuple (it[0], it[1][1], it[2][1], it[3][1], it[4][1])
              }     
       .set { to_mix }
separated_by_segment_toBeImputed_and_toBeUsedAsImputationPanel =
       to_mix.ch_one.mix(to_mix.ch_two)

process phasing {
 publishDir params.phasedDir, mode: 'move'

 //cpus 8 //8 necessary, but optimal value is 2
 cpus params.cpus
 label 'Phasing'
 tag "${intervalname}.${input}"

 input:
 tuple val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment_toBeImputed_and_toBeUsedAsImputationPanel
 
 output:
 set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.phased.vcf.gz"), file("${remExt(vcf.name)}.phased.vcf.gz.tbi") into separated_by_segment_toBeImputed_and_toBeUsedAsImputationPanel_phased

 script:
 chr = intervalname.split('_')[0]
 """
 uname -a | awk '{print \$2}'
 # Phasing
 ${params.refDir}/Eagle_v2.4.1/eagle \
          --vcf ${vcf} \
          --chrom  ${chr} \
          --geneticMapFile ${params.refDir}/imputation/mapChr/eagle_${chr}_b38.map \
          --numThreads=${params.cpus} \
          --Kpbwt=20000 \
          --outPrefix ${remExt(vcf.name)}.phased
 bcftools index -t ${remExt(vcf.name)}.phased.vcf.gz
 """
}

// Separate phased channel into one for bref
toBeImputed = Channel.create()
imputationPanel = Channel.create()
separated_by_segment_toBeImputed_and_toBeUsedAsImputationPanel_phased
       .choice(toBeImputed, imputationPanel) { it[2] == remPath(params.toBeImputed) ? 0 : 1 }

process bref_imp_panel {
       label 'bref'
       tag "${intervalname}.${input}"

       input:
       tuple val(order), val(intervalname), val(input), file(vcf), file(idx) from imputationPanel
       
       output:
       set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.bref"), file("equaliser_element") into imputationPanel_bref
       
       script:
       """
       java -jar ${params.refDir}/bref.27Jan18.7e1.jar ${vcf}
       touch equaliser_element
       """
}

// Combine toBeImputed and ImputationPanel channels
imputation_ch = toBeImputed.map {tuple (it,0)}
       .mix(imputationPanel_bref.map {tuple (it,1)})
       .toSortedList({ a,b -> a[1] <=> b[1] })
       .flatten().buffer (size: 6)
       .map { tuple(it[0..4]) }
       .groupTuple(by:[0,1])

// Customise manipulation steps
process manipulate_segment_imputation {
 //publishDir = params.publishDir
 cpus params.cpus

 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx) from imputation_ch

 output:
 set val(order), val(intervalname), val(input), file("${output}.INFO.vcf.gz"), file("${output}.INFO.vcf.gz.tbi") into segments_ready_for_collection_imputed

 script:
 chr = intervalname.split('_')[0]
 start = intervalname.split('_')[1]
 stop = intervalname.split('_')[2]
 output = "${input[0]}.imputed_with.${input[1]}.${intervalname}"
 input = input[0] // for easier reproducibility
 cp = 30
 """
 # Imputation
 java -Xss5m -Xmx64g -jar ${params.refDir}/beagle.27Jan18.7e1.jar \
          gt=${vcf[0]} \
          ref=${vcf[1]} \
          map=${params.refDir}/imputation/Imputation/dockers/reference-data-full/reference-data/map/beagle_${chr}_b38.map \
          out=${output} \
          chrom=${chr}:${start}-${stop} \
          window=500000 \
          nthreads=${params.cpus} \
          niterations=10 \
          gprobs=true \
          ne=20000 \
          impute=true \
          seed=-99999
 # Add impute2 like INFO score
 bcftools index -t ${output}.vcf.gz
 bcftools +fill-tags ${output}.vcf.gz -- -t AF,AC |
 bcftools +impute-info -Oz -o ${output}.INFO.vcf.gz
 bcftools index -t ${output}.INFO.vcf.gz
 """
}

segments_sample_ready_for_collection_collected = segments_ready_for_collection_imputed
 .toSortedList({ a,b -> a[0] <=> b[0] })
 .flatten().buffer ( size: 5 )
 .groupTuple(by:[2])
 
// Concatanate segments
process concatanate_segments {
 publishDir params.publishDir, mode: 'move', overwrite: true

 input:
 set val(order), val(intervalname), val(input), file(vcf_all), file(idx_all) from segments_sample_ready_for_collection_collected 
 output:
 file output
 script:
 output = "${vcf_all[0].name}" - "${intervalname[0]}."
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 bcftools concat --naive -f vcfFiles.txt -Oz -o ${output}
 
 """
}
/*
*/