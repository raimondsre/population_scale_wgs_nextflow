#!/usr/bin/env nextflow
//nextflow.enable.dsl = 2

params.input = 'Hello world!'
process convertToUpper {
       publishDir '.'
       output:
       file "test"
       """
       echo ${params.input} > test
       """
}