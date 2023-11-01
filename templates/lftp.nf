#!/usr/bin/env nextflow

//before executing, set up lftp authentification

params.sampleLocation = './for.trimming'
params.batchDir = '.'
//params.batchName = 'lv_reference_20220722_502samples'
//params.trimGaloreContainer = '/mnt/beegfs2/beegfs_large/raimondsre_add2/genome_analysis/trim_galore_0.6.7.sif'
params.hpc_billing_account = 'bmc_flpp_0151'

// Read the input file containing sample ID and it's location
Channel
       .fromPath(params.sampleLocation)
       .splitCsv(header:false, sep:'\t',strip:true)
       .map { row -> tuple(row[0], row[1]) }
       //.filter({it[1].contains('25408')})
       .set { for_lftp }
def remPath(String fileName) {return fileName.replaceAll(/.*\//,'')}

//for_lftp.subscribe{ println it }

process file_transfer { 
       publishDir params.batchDir, mode: 'move', overwrite: false
       cpus 1
       //executor 'pbs'
       clusterOptions "-l nodes=wn61 -A ${hpc_billing_account}"

       input:
       set val(read1_lftp), val(read2_lftp) from for_lftp
       
       output:
       set file(read1), file(read2)

       script:
       read1 = remPath(read1_lftp)
       read2 = remPath(read2_lftp)
       """
       lftp -e "set ssl:verify-certificate no; set net:connection-limit 2; get ${read1_lftp} -o ${read1} & get ${read2_lftp} -o ${read2} & wait; exit"
       wait
       """
}

/*
*/