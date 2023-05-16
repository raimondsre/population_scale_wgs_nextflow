#!/usr/bin/env nextflow
//nextflow.enable.dsl = 2

params.str = 'Hello world!'
process splitLetters {
       output:
       file 'chunk_*' into letters mode flatten
       """
       printf '${params.str}' | split -b 6 - chunk_
       """
}
letters
       .into{ letters; letters2 }

process convertToUpper {
       input:
       file x from letters 
       output:
       file ("file") into result
       """
       echo $x > file
       """
}

process testing_template {
       publishDir '.', mode: 'copy', overwrite: true
       input:
       file input from letters2
       output:
       file "test" 

       shell:
       '''
       nextflow run !{projectDir}/test_simple_template.nf --input ${PWD}/!{input}
       '''
}

result.subscribe {
       println it
}