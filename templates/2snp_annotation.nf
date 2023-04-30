#!/usr/bin/env nextflow
// main_bySegmentOnly.nf
params.publishDir = './results'

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
 .filter({it[1].contains('chr1')})
 .filter({it[2].contains('240000001')})
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
       # bcftools view ${vcf} ${chr}:${start}-${stop} |
       # bcftools view --exclude 'POS<${start}' |
       # bcftools view --exclude 'POS>${stop}' -Oz -o ${input}.${intervalname}.vcf.gz
       # bcftools index -t ${input}.${intervalname}.vcf.gz
       touch ${input}.${intervalname}.vcf.gz
       touch ${input}.${intervalname}.vcf.gz.tbi
       variantsPresent=1
       # if [ `bcftools view ${input}.${intervalname}.vcf.gz --no-header | wc -l` -eq 0 ]; then variantsPresent=0; fi
 """
}
separated_by_segment = separated_by_segment.filter { it[5] == "1"  }.map {it - it[5]}

separated_by_segment.subscribe {println it}


// Customise manipulation steps
process manipulate_segment_vep {
 publishDir params.publishDir
 //cpus 1

 input:
 set val(order), val(intervalname), val(input), file(vcf), file(idx) from separated_by_segment

 output:
 set val(order), val(intervalname), val(input), file("${intervalname}.all.vep.counted"), file("${remExt(vcf.name)}.vep") into segments_ready_for_collection

 script:
 vcf_name = vcf.name
 """
 mkdir vcf_file
 cp ${vcf} vcf_file/
 singularity run /home_beegfs/raimondsre/programmas/vep.sif vep --offline \
    --dir_cache /home/raimondsre/.vep --species homo_sapiens --vcf --assembly GRCh38 \
    --af_gnomade --variant_class --biotype --check_existing --symbol --compress_output bgzip \
    --custom /home/raimondsre/.vep/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNDN \
    --canonical \
    -i vcf_file/${vcf_name} \
    -o ${remExt(vcf.name)}.vep.vcf.gz
 # VCF to txt
 # If annotated column is added (e.g. PHENO), you have to add corresponding column to countVEPfeatures.R file line 21!
 bcftools +split-vep -d -f '%ID %VARIANT_CLASS %CLIN_SIG %Consequence %Existing_variation %gnomADe_AF %ClinVar_CLNDN\n' \
    ${remExt(vcf.name)}.vep.vcf.gz > ${remExt(vcf.name)}.vep
 # Count features
 Rscript ${projectDir}/countVEPfeatures.R --input ${remExt(vcf.name)}.vep --interval ${intervalname} --sample all --original_file_name ${input}
 """
}

// Arrange segments and group by input file name
segments_ready_for_collection_collected = segments_ready_for_collection
 //sorts by the first variable [0]
 .toSortedList({ a,b -> a[0] <=> b[0] })
 //flattens whole channel and then remakes sets by 4 elements (if input sets were 5 elements, number should be changed)
 .flatten().buffer ( size: 4 )
 //squishes all elements into one set, collects. For each set variable [2] there will be one channel element.
 .groupTuple(by:[2])

// Concatanate segments
process concatanate_segments {
 publishDir params.publishDir, mode: 'move', overwrite: true
 //cpus 16
 input:
 set val(order), val(intervalname), val(input), file(vep_all) from segments_ready_for_collection_collected 
 output:
 set file("${input}.vep.counted")
 script:
 """
 cat ${vep_all.join(' ')} > ${input}.vep.counted
 """
}
*/