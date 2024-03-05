#!/usr/bin/env nextflow
//    nextflow run raimondsre/population_scale_wgs_nextflow/test_simple.nf -r main -latest --input /home/raimondsre/test/input.csv --output /home/raimondsre/test/output.csv
nextflow.enable.dsl=2

params.input = './input.csv'
params.output = './output.csv'

process calculate_polygenic_score {
       executor = 'local'
       input:
       path input_file
       """
       cat ${input_file} > ${params.output}
       """
}

workflow {
       input_file = file(params.input)

       input_file | calculate_polygenic_score
}