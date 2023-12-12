#!/usr/bin/env nextflow
// Workflow to performs 
//     1. Input 23andMe or VCF format genome harmonisation to GRCh38
//     2. Quality control
//     3. Imputation
//     4. PGS calculation
//     5. Visualisation and/or JSON preparation
// Workflow to impute VCF with another VCF (reference panel)
// Outputs imputed VCF and variants counted by AF and INFO score

params.publishDir = './results'
params.refDir = '/home_beegfs/groups/bmc/genome_analysis_tmp/hs/ref'
params.phasedDir = '/mnt/beegfs2/home/groups/bmc/references/populationVCF/phased' // Contains phased and bref corrected segments
params.hg37fasta = '/home_beegfs/groups/bmc/genome_analysis_tmp/hs/ref/human_g1k_v37.fasta'
params.hg38fasta = '/home_beegfs/groups/bmc/genome_analysis_tmp/hs/ref/Homo_sapiens_assembly38.fasta'
params.cpus = 8

params.toBeImputed = './'
params.imputationPanel1 = './'
//params.VCFfile = './merged.two.vcf.gz'
//params.intervalsBed = './hg38chr25int5e6.bed'

// Define channels for intervals and initial .vcf.gz file
// Input file to be imputed
Channel
 .fromPath(params.toBeImputed)
 .map { file(it) }
 .into { input_genome; vcf_toBeImputed_extractSamples }
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
 //.filter({!it[4].contains('chr5_10000001_15000000')}) // Imputation problematic with the following segments: chr9_40000001_45000000,
 .into { intervals1; intervals2 }
 
// Samples in input VCF
// process extract_vcf_samples {
//  input:
//  tuple file(vcf), file(idx) from vcf_toBeImputed_extractSamples
//  output:
//  file 'samples' into samples_ch mode flatten
//  script:
//  """
//  bcftools query -l ${vcf} > samples
//  """
// }
// counter2 = 0
// samples_ch
//  .splitText() {it.replaceFirst(/\n/,'')}
//  .map {value ->
//         counter2 += 1
//         [counter2, value].flatten()}
// .into { samples_ch1; samples_ch2 }
// Define function to remove .vcf.gz extension
def remPath(String fileName) {return fileName.replaceAll(/.*\//,'').replaceFirst(/\.vcf\.gz$/,'')}
def remExt(String fileName) {return fileName.replaceFirst(/\.vcf\.gz$/,'')}
def remExtBref(String fileName) {return fileName.replaceFirst(/\.bref$/,'')}
def getExt(String fileName) {return fileName.replaceAll(/.*\./,'')}

//###
//### Analysis
//###

// Harmonise genomes
process harmonisation {
       publishDir params.publishDir, mode: 'copy', overwrite: true
       
       input:
       file(genome) from input_genome
       output:
       set file("merged.vcf.gz")

       script:
       intput_ext = getExt(genome.name)
       """
       if [ ${intput_ext} == "zip" ]; then
              unzip ${genome}
              mv *txt genome.txt
              plink --23file genome.txt --keep-allele-order --output-chr MT --snps-only just-acgt --recode vcf --out genome
              #normalise to hg19
              bcftools +fixref genome.vcf -Oz -o genome.vcf.gz -- -f ${params.hg37fasta} -m top
              plink --vcf genome.vcf.gz --keep-allele-order --output-chr chrM --recode vcf --out genome.chrM
              export _JAVA_OPTIONS="-Xmx16g"
              ~/.conda/envs/picard/bin/picard LiftoverVcf \
                     I=genome.chrM.vcf \
                     O=output.hg38.vcf \
                     CHAIN=/home/raimondsre/array/input_data/ref/hg19ToHg38.over.chain.gz \
                     REJECT=rejected_variants.vcf \
                     R=/home_beegfs/groups/bmc/genome_analysis_tmp/hs/ref/Homo_sapiens_assembly38.fasta
              bcftools annotate --set-id +'%CHROM:%POS:%REF:%ALT' output.hg38.vcf -Oz -o output.hg38.vcf.gz
              bcftools +fixref output.hg38.vcf.gz -Oz -o output.hg38.altFilter.fixref.vcf.gz -- -f ${params.hg38fasta} -m top
              bcftools index -t output.hg38.altFilter.fixref.vcf.gz
              bcftools merge output.hg38.altFilter.fixref.vcf.gz /home_beegfs/groups/bmc/genome_analysis_tmp/hs/analysis/pgr_kalkulators/nextflow/gsa.array.192.af_filt.fixref.vcf.gz -Oz -o merged.vcf.gz
       fi
       touch normalised_genome.vcf.gz
       """
}

/*

// Make single channel for intervals and vcf file to be imputed
vcfIntervals_toBeImputed = intervals1.combine(vcf_toBeImputed)
// Make single channel for intervals and vcf file to be used as imputation panel
vcfIntervals_toBeUsedAsImputationPanel = intervals2.combine(vcf_imputation_panel)
// Concatanate to-be-separated channels into a single channel
vcfIntervals_toBeImputed_and_toBeUsedAsImputationPanel = vcfIntervals_toBeImputed.mix(vcfIntervals_toBeUsedAsImputationPanel)

// Separate VCF into fragments (has to be before separating by sample)
process separateVCF {
       publishDir params.phasedDir, mode: 'copy', overwrite: false
       //publishDir = params.publishDir
       input:
       tuple val(order), val(chr), val(start), val(stop), val(intervalname), file(vcf), file(idx) from vcfIntervals_toBeImputed_and_toBeUsedAsImputationPanel
       
       output:
       set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi"), env(variantsPresent) into separated_by_segment_toBeImputed_and_toBeUsedAsImputationPanel

       script:
       input = remExt(vcf.name)  
       """
       if [ -e ${params.phasedDir}/${input}.${intervalname}.vcf.gz ]; then
              cp ${params.phasedDir}/${input}.${intervalname}.vcf.gz ${input}.${intervalname}.vcf.gz
              cp ${params.phasedDir}/${input}.${intervalname}.vcf.gz.tbi ${input}.${intervalname}.vcf.gz.tbi
       else
              bcftools view ${vcf} ${chr}:${start}-${stop} |
              bcftools view --exclude 'POS<${start}' |
              bcftools view --exclude 'POS>${stop}' |
              bcftools norm --multiallelics - |
              bcftools norm --remove-duplicates |
              bcftools view -Oz -o ${input}.${intervalname}.vcf.gz
              bcftools index -t ${input}.${intervalname}.vcf.gz
       fi
       variantsPresent=1
       # Check wether VCF segment has variants present
       if [ `bcftools view ${input}.${intervalname}.vcf.gz --no-header | wc -l` -lt 50 ]; then variantsPresent=0; fi
       # Check wether genetic map file contains variants in segment
       if [ `awk '\$2 >= ${start} && \$2 <= ${stop} {print \$0}' ${params.refDir}/imputation/mapChr/eagle_${chr}_b38.map | wc -l` -lt 50 ]; then variantsPresent=0; fi
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
 publishDir params.phasedDir, mode: 'copy', overwrite: false

 //cpus 8 //8 necessary, but optimal value is 2
 cpus params.cpus
 label 'Phasing'
 tag "${intervalname}.${input}"

 input:
 tuple val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment_toBeImputed_and_toBeUsedAsImputationPanel
 
 output:
 set val(order), val(intervalname), val(input), file("${phased}.vcf.gz"), file("${phased}.vcf.gz.tbi") into separated_by_segment_toBeImputed_and_toBeUsedAsImputationPanel_phased

 script:
 chr = intervalname.split('_')[0]
 phased = remExt(vcf.name)+".phased"
 """
 uname -a | awk '{print \$2}'
 # Phasing
 # If phased segment already exists, take it. Otherwise phase anew
 if [ -e ${params.phasedDir}/${phased}.vcf.gz ]; then
  cp ${params.phasedDir}/${phased}.vcf.gz ${phased}.vcf.gz
  cp ${params.phasedDir}/${phased}.vcf.gz.tbi ${phased}.vcf.gz.tbi
 else
  ${params.refDir}/Eagle_v2.4.1/eagle \
          --vcf ${vcf} \
          --chrom  ${chr} \
          --geneticMapFile ${params.refDir}/imputation/mapChr/eagle_${chr}_b38.map \
          --numThreads=${params.cpus} \
          --Kpbwt=20000 \
          --outPrefix ${phased}
  bcftools index -t ${phased}.vcf.gz
 fi
 """
}

// Separate phased channel into one for bref
toBeImputed = Channel.create()
imputationPanel = Channel.create()
separated_by_segment_toBeImputed_and_toBeUsedAsImputationPanel_phased
       .choice(toBeImputed, imputationPanel) { it[2] == remPath(params.toBeImputed) ? 0 : 1 }

process bref_imp_panel {
       publishDir params.phasedDir, mode: 'copy', overwrite: false
       label 'bref'
       tag "${intervalname}.${input}"

       input:
       tuple val(order), val(intervalname), val(input), file(vcf), file(idx) from imputationPanel
       
       output:
       set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.bref"), file("equaliser_element") into imputationPanel_bref
       
       script:
       """
       if [ -e ${params.phasedDir}/${remExt(vcf.name)}.bref ]; then
        cp ${params.phasedDir}/${remExt(vcf.name)}.bref ${remExt(vcf.name)}.bref
        touch equaliser_element
       else
        java -jar ${params.refDir}/bref.27Jan18.7e1.jar ${vcf}
        touch equaliser_element
       fi
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
(segments_ready_for_collection_imputed, segments_ready_for_collection_imputed_for_info_counting) = segments_ready_for_collection_imputed.into(2)
segments_sample_ready_for_collection_collected = segments_ready_for_collection_imputed
 .toSortedList({ a,b -> a[0] <=> b[0] })
 .flatten().buffer ( size: 5 )
 .groupTuple(by:[2])
 
// Concatanate segments
process concatanate_segments {
 publishDir params.publishDir, mode: 'copy', overwrite: true

 input:
 set val(order), val(intervalname), val(input), file(vcf_all), file(idx_all) from segments_sample_ready_for_collection_collected 
 output:
 set file(output), file("${output}.tbi")
 script:
 output = "${vcf_all[0].name}" - "${intervalname[0]}."
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 bcftools concat -f vcfFiles.txt -Oz -o ${output}
 bcftools index -t ${output}
 """
}
/*
// Counting variant number by info score
process count_by_info_score {
       //publishDir params.publishDir, mode: 'copy', overwrite: true
       input:
       set val(order), val(intervalname), val(input), file(vcf), file(idx) from segments_ready_for_collection_imputed_for_info_counting

       output:
       set val(name), file(output_counted) into counted_segments_ready_for_collection

       script:
       output_full = remExt(vcf.name)+".txt"
       output_counted = remExt(vcf.name)+".counted.txt"
       name = "${remExt(vcf.name)}" - ".${intervalname}.INFO"
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