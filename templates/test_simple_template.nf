#!/usr/bin/env nextflow
//nextflow.enable.dsl = 2

params.input = 'Hello world!'
process convertToUpper {
       output:
       file "test"
       """
       echo ${params.input} > test
       """
}