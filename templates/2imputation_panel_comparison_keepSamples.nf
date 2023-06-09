#!/usr/bin/env nextflow
// Pipeline to impute VCF with another VCF (reference panel)
// Outputs imputed VCF and variants counted by AF and INFO score

// Imputation pipeline to variably select subsets of toBeImputed set
// Here is an example script to launch imputation for multiple subsets of VCF file:
// samples being two column, nominatind sample ID and groupl (e.g. country)

//                 samples=/home_beegfs/groups/bmc/references/phenotypes/aadr_european_10plus.samples
//                 samplesDir=/home_beegfs/groups/bmc/genome_analysis_tmp/hs/analysis/lv_genome_reference_20220722/imputation/imputing_different_european_populations
                
//                 cut -f2 $samples | sort | uniq | cat | while read in; do
//                 country=$in
//                 awk -v var=$country '$2 == var {print $1}' $samples > $samplesDir/samples.$country
                
//                 set=v54.1.p1_HO_public.hg38.normalised.${country}
//                 toImpute=/home_beegfs/groups/bmc/references/populationVCF/original/aadr/v54.1.p1_HO_public.hg38.normalised.vcf.gz
//                 panel=/home_beegfs/groups/bmc/genome_analysis_tmp/hs/analysis/imp/20220722/NYpanel.vcf.gz
//                 #panel_id=lvbmc_502
//                 panel_id=1000G_2022
                
//                 refDir=/mnt/beegfs2/home/groups/bmc/genome_analysis_tmp/hs/ref
//                 intervals=${refDir}/intervals5mil
//                 config=${refDir}/configurations/nextflow_imputation_noClean.config
//                 baseDir=/home_beegfs/groups/bmc/genome_analysis_tmp/hs/analysis/lv_genome_reference_20220722
//                 resultsDir=${baseDir}/imputation/imputing_different_european_populations
//                 runID=${set}_set_${panel_id}_panel
//                 workDir=${resultsDir}/${runID}; mkdir -p $workDir; cd $workDir
//                 impPanel=${panel}
//                 nextflow run raimondsre/population_scale_WGS_nextflow/templates/2imputation_panel_comparison_keepSamples.nf \
//                 -r main -latest \
//                 --toBeImputed ${toImpute} \
//                 --imputationPanel1 ${impPanel} \
//                 --publishDir ${resultsDir} \
//                 --intervalsBed $intervals \
//                 -c $config \
//                 --subsetname ${country} \
//                 --samplesToKeep $samplesDir/samples.${country}
                
//                 done

params.publishDir = './results'
params.refDir = '/home_beegfs/groups/bmc/genome_analysis_tmp/hs/ref'
params.phasedDir = '/mnt/beegfs2/home/groups/bmc/references/populationVCF/phased' // Contains phased and bref corrected segments
params.cpus = 4
params.maxforks = 32
params.executor = 'local'

params.toBeImputed = './'
params.imputationPanel1 = './'
params.samplesToKeep = './keep.samples'
params.subsetname = 'no_subset'

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
 //.filter({it[4].contains('chr4_190000001_190214555')}) // Imputation problematic with the following segments: chr9_40000001_45000000,
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
//vcfIntervals_toBeImputed_and_toBeUsedAsImputationPanel = vcfIntervals_toBeImputed.mix(vcfIntervals_toBeUsedAsImputationPanel)

//###
//### Analysis
//###

// Separate VCF into fragments (has to be before separating by sample)
process separateVCF_imputation_panel {
       executor params.executor
       publishDir params.phasedDir, mode: 'symlink', overwrite: false
       //publishDir = params.publishDir
       input:
       tuple val(order), val(chr), val(start), val(stop), val(intervalname), file(vcf), file(idx) from vcfIntervals_toBeUsedAsImputationPanel
       
       output:
       set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi"), env(variantsPresent) into separated_by_segment_toBeUsedAsImputationPanel

       script:
       input = remExt(vcf.name) 
       """
       if [ -e ${params.phasedDir}/${input}.${intervalname}.vcf.gz ]; then
              ln -s ${params.phasedDir}/${input}.${intervalname}.vcf.gz ${input}.${intervalname}.vcf.gz
              ln -s ${params.phasedDir}/${input}.${intervalname}.vcf.gz.tbi ${input}.${intervalname}.vcf.gz.tbi
       else
              bcftools view ${vcf} ${chr}:${start}-${stop} |
              bcftools view --exclude 'POS<${start}' |
              bcftools view --exclude 'POS>${stop}' |
              bcftools norm --remove-duplicates |
              bcftools view -c3 -Oz -o ${input}.${intervalname}.vcf.gz
              bcftools index -t ${input}.${intervalname}.vcf.gz
       fi
       variantsPresent=1
       # Check wether VCF segment has variants present
       if [ `bcftools view ${input}.${intervalname}.vcf.gz --no-header | wc -l` -lt 50 ]; then variantsPresent=0; fi
       # Check wether genetic map file contains variants in segment
       if [ `awk '\$2 >= ${start} && \$2 <= ${stop} {print \$0}' ${params.refDir}/imputation/mapChr/eagle_${chr}_b38.map | wc -l` -eq 0 ]; then variantsPresent=0; fi
       """
}
process separateVCF_to_be_imputed {
       executor params.executor
       //publishDir params.phasedDir, mode: 'copy', overwrite: false
       //publishDir = params.publishDir
       input:
       tuple val(order), val(chr), val(start), val(stop), val(intervalname), file(vcf), file(idx) from vcfIntervals_toBeImputed
       
       output:
       set val(order), val(intervalname), val(input), file("${input}.${intervalname}.vcf.gz"), file("${input}.${intervalname}.vcf.gz.tbi"), env(variantsPresent) into separated_by_segment_toBeImputed

       script:
       input = remExt(vcf.name)+"."+params.subsetname
       """
       bcftools view ${vcf} ${chr}:${start}-${stop} |
              bcftools view --exclude 'POS<${start}' |
              bcftools view --exclude 'POS>${stop}' |
              bcftools norm --remove-duplicates |
              bcftools view -c3 -Oz -o ${input}.${intervalname}.vcf.gz
       if [ ${params.subsetname} != "no_subset" ]; then
              bcftools view -S ${params.samplesToKeep} --force-samples ${input}.${intervalname}.vcf.gz -Oz -o temp.vcf.gz
              mv temp.vcf.gz ${input}.${intervalname}.vcf.gz
       fi
       bcftools index -t ${input}.${intervalname}.vcf.gz

       variantsPresent=1
       # Check wether VCF segment has variants present
       if [ `bcftools view ${input}.${intervalname}.vcf.gz --no-header | wc -l` -lt 50 ]; then variantsPresent=0; fi
       # Check wether genetic map file contains variants in segment
       if [ `awk '\$2 >= ${start} && \$2 <= ${stop} {print \$0}' ${params.refDir}/imputation/mapChr/eagle_${chr}_b38.map | wc -l` -eq 0 ]; then variantsPresent=0; fi
       """
}


// Select toBeImputed - impPanel segment pairs only, if both have at least 1 variant after filtering
separated_by_segment_toBeUsedAsImputationPanel
       .mix(separated_by_segment_toBeImputed)
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
executor params.executor

 publishDir params.phasedDir, mode: 'symlink', overwrite: false
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
  ln -s ${params.phasedDir}/${phased}.vcf.gz ${phased}.vcf.gz
  ln -s ${params.phasedDir}/${phased}.vcf.gz.tbi ${phased}.vcf.gz.tbi
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
       .choice(toBeImputed, imputationPanel) { it[2] == remPath(params.toBeImputed)+"."+params.subsetname  ? 0 : 1 }

process bref_imp_panel {
       executor params.executor
       publishDir params.phasedDir, mode: 'symlink', overwrite: false
       label 'bref'
       tag "${intervalname}.${input}"

       input:
       tuple val(order), val(intervalname), val(input), file(vcf), file(idx) from imputationPanel
       
       output:
       set val(order), val(intervalname), val(input), file("${remExt(vcf.name)}.bref"), file("equaliser_element") into imputationPanel_bref
       
       script:
       """
       if [ -e ${params.phasedDir}/${remExt(vcf.name)}.bref ]; then
        ln -s ${params.phasedDir}/${remExt(vcf.name)}.bref ${remExt(vcf.name)}.bref
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
 bcftools view -c1 |
 bcftools +impute-info -Oz -o ${output}.INFO.vcf.gz
 bcftools index -t ${output}.INFO.vcf.gz
 """
}
(segments_ready_for_collection_imputed, segments_ready_for_collection_imputed_for_info_counting) = segments_ready_for_collection_imputed.into(2)
/*
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
 file output 
 script:
 output = "${vcf_all[0].name}" - "${intervalname[0]}."
 """
 echo "${vcf_all.join('\n')}" > vcfFiles.txt
 bcftools concat -f vcfFiles.txt -Oz -o ${output}
 
 """
}
*/
// Counting variant number by info score
process count_by_info_score {
       executor params.executor
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
       executor params.executor
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