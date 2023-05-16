#!/usr/bin/env nextflow
//nextflow.enable.dsl = 2

params.input = 'Hello world!'
process convertToUpper {
       output:
       stdout result
       """
       echo ${params.input}
       """
}
result.subscribe {
       println it.trim()
}