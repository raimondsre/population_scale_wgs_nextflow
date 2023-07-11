#!/usr/bin/env nextflow

params.input = './input.csv'
params.output = './output.csv'

process calculate_polygenic_score {
       """
       cat ${params.input} > ${params.output}
       """
}