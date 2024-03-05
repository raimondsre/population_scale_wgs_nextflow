#!/usr/bin/env nextflow

params.input = './input.csv'
params.output = './output.csv'

process calculate_polygenic_score {
       executor = 'local'
       """
       cat ${params.input} > ${params.output}
       Exit 1
       """
}