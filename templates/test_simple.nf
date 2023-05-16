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
       stdout result
       """
       echo $x 
       """
}

process testing_template {
       input:
       file x from letters2
       output:
       stdout result2

       script:
       template 'test_simple_template.nf --input ${x}'
}

result.subscribe {
       println it.trim()
}