#!/usr/bin/env nextflow

params.sampleLocation = './for.trimming'
params.batchDir = '.'
params.batchName = 'lv_reference_20220722_502samples'
params.EGAencryptor = '/home/raimondsre/programms/EGA-Cryptor-2.0.0/ega-cryptor-2.0.0.jar'

// Read the input file containing sample ID and it's location
Channel
       .fromPath(params.sampleLocation)
       .splitCsv(header:false, sep:'\t',strip:true)
       .map { row -> tuple(row[0], row[1], row[2]) }
       .multiMap { 
              ch_one: tuple (it[0], it[1], 1)
              ch_two: tuple (it[0], it[2], 2)
              }
       .set { to_mix }
to_mix.ch_one.mix(to_mix.ch_two)
       .into { to_encrypt; intervals2 }
//to_encrypt.subscribe { println it}

// Encrypt each file
process ega_encrypt {
       publishDir params.batchDir+"/"+params.batchName, mode: 'move', overwrite: false
       cpus 4
       
       input:
       set val(SAMPLE_ID), path(read), val(read_num) from to_encrypt

       output:
       set file(read_encrypted), file(read_encrypted_checksum), file(read_unencrypted_checksum) into read_encrypted

       // when:
       // (params.batchDir+"/"+params.batchName+"/"+read.name()+".gpg").isEmpty()

       script:
       read_encrypted = read+".gpg"
       read_encrypted_checksum = read+".gpg.md5"
       read_unencrypted_checksum = read+".md5"
       """
       if [ ! -f ${params.batchDir}/${params.batchName}/${read_encrypted} ]; then
       java -jar ${params.EGAencryptor} -i ${read}
       mv output-files/* .
       else
       ln -s ${params.batchDir}/${params.batchName}/${read_encrypted} > ${read_encrypted}
       ln -s ${params.batchDir}/${params.batchName}/${read_encrypted_checksum} > ${read_encrypted_checksum}
       ln -s ${params.batchDir}/${params.batchName}/${read_unencrypted_checksum} > ${read_unencrypted_checksum}
       fi
       """
}

// Upload files to EGA
process ega_upload {
       cpus 1
       
       input:
       set val(read_encrypted), path(read_encrypted_checksum), val(read_unencrypted_checksum) from read_encrypted

       script:
       """
       lftp ftp.ega.ebi.ac.uk -e "put ${read_encrypted} -o lv_reference_20220722_502samples/; put ${read_encrypted_checksum} -o lv_reference_20220722_502samples/; put ${read_unencrypted_checksum} -o lv_reference_20220722_502samples/; exit"
       """
}

/*
*/