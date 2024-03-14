#!/usr/bin/env nextflow
//    nextflow run raimondsre/population_scale_wgs_nextflow/test_simple.nf -r main -latest --input /home/raimondsre/test/input.csv --output /home/raimondsre/test/output.csv

params.input = './input.csv'
params.output = './output.csv'

process calculate_polygenic_score {
       executor = 'local'
       container '/home/raimondsre/analysis/hs/genome/prs/app/continental_ethnicity/picard:3.1.1--hdfd78af_0'
       input:
       path input_file
       """
       if grep -q error "${input_file}"; then
              echo "Test cases failed"
              exit 1
       fi
       cat ${input_file} > ${params.output}
       """
}

workflow {
       Channel.fromPath(params.input).set{ test_file }

       test_file | calculate_polygenic_score
}