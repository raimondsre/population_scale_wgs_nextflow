#!/usr/bin/env nextflow

params.input = './input.csv'

process calculate_polygenic_score {
       publishDir "./results"
       output:
       file "output.csv" 
       """
       cat ${params.input} > output.csv
       """
