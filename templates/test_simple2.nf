#!/usr/bin/env nextflow
//    nextflow run raimondsre/population_scale_wgs_nextflow/templates/test_simple2.nf -r main -latest --input /home/raimondsre/test/input.csv --output /home/raimondsre/test/output.csv

params.input = './input.csv'
params.output = './output.csv'

process calculate_polygenic_score {
       executor = 'local'
       input:
       path input_file
       """
       if grep -q "Pipeline completed with errors" "${input_file}"; then
              echo "Test cases failed"
              exit 1
       else
              echo "Test 2 case did not fail"
       fi
       cat ${input_file} >> ${params.output}
       """
}

workflow {
       Channel.fromPath(params.input).set{ test_file }

       test_file | calculate_polygenic_score
}