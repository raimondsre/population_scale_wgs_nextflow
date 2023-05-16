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
       publishDir '.'
       input:
       file ("file") from letters2
       output:
       file "test" 

       shell:
       '''
       nextflow run !{projectDir}/test_simple_template.nf --input file
       '''
}

result.subscribe {
       println it
}