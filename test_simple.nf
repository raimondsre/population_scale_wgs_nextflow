#!/usr/bin/env nextflow
//    nextflow run raimondsre/population_scale_wgs_nextflow/test_simple.nf -r main -latest --input /home/raimondsre/test/input.csv --output /home/raimondsre/test/output.csv

params.input = './input.csv'
params.output = './output.csv'

process calculate_polygenic_score {
       executor = 'pbs'
       input:
       path input_file
       """
       cat ${input_file} > ${params.output}
       """
}

workflow {
       Channel.fromPath(params.input).set{ test_file }

       test_file | calculate_polygenic_score
}