#!/usr/bin/env nextflow
//    nextflow run raimondsre/population_scale_wgs_nextflow/templates/test_simple.nf -r main -latest --input /home/raimondsre/test/input.csv --output /home/raimondsre/test/output.csv
nextflow.enable.dsl=1

params.input = './input.csv'
params.output = './output.csv'

process calculate_polygenic_score {
       executor = 'local'
       """
       cat ${params.input} > ${params.output}
       exit 1
       """
}