#!/usr/bin/env nextflow
//    nextflow run raimondsre/population_scale_wgs_nextflow/test_simple.nf -r main -latest --input /home/raimondsre/test/input.csv --output /home/raimondsre/test/output.csv

params.input = './input.csv'
params.output = './output.csv'

process calculate_polygenic_score {
       executor = 'local'
       input:
       path input_file
       """
       if grep -q error "${input_file}"; then
              echo "Test cases failed"
              exit 1
       fi
       cat ${input_file} > ${params.output}
       
       nextflow run raimondsre/population_scale_wgs_nextflow/templates/test_simple2.nf -r main -latest \
       --input /home_beegfs/groups/bmc/genome_analysis_tmp/hs/analysis/pgr_kalkulators/nextflow/input.csv \
       --output /home_beegfs/groups/bmc/genome_analysis_tmp/hs/analysis/pgr_kalkulators/nextflow/output.csv \
       -c /home_beegfs/groups/bmc/genome_analysis_tmp/hs/analysis/pgr_kalkulators/nextflow/nextflow_imputation.config
       """
}

workflow {
       Channel.fromPath(params.input).set{ test_file }

       test_file | calculate_polygenic_score
}