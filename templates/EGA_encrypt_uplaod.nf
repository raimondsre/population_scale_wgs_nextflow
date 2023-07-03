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
       //.filter({it[1].contains('25408')})
       .multiMap { 
              ch_one: tuple (it[0], it[1], 1)
              ch_two: tuple (it[0], it[2], 2)
              }
       .set { to_mix }
to_mix.ch_one.mix(to_mix.ch_two)
       // Filter out fastq files that does not exist
       .filter { SAMPLE_ID, read, read_num ->
        def fastq_filename = read.tokenize('/').last()
        def file = new File(read)
        def file_exists = file.exists()
        if (! file_exists) println ">>> WARNING: File ${read} does not exist and will not be included"
        file_exists
       }
       // Filter out fastq files that has already been uploaded
       .filter { SAMPLE_ID, read, read_num ->
        def fastq_path = read
        def fastq_filename = fastq_path.tokenize('/').last()
        def file = new File("${params.batchDir}/${params.batchName}/list_of_ega_uploaded/${fastq_filename}.gpg")
        def uploaded_to_ega = file.exists()
        if (uploaded_to_ega) prinltn ">>> WARNING: FILE ${fastq_filename} of ${SAMPLE_ID} sample already encrypted and uploaded"
        !uploaded_to_ega
       }
       .into { to_encrypt; intervals2 }

// Encrypt each file
process ega_encrypt {
       publishDir params.batchDir+"/"+params.batchName, mode: 'move', overwrite: false
       cpus 4
       
       input:
       set val(SAMPLE_ID), path(read), val(read_num) from to_encrypt

       output:
       set file(read_encrypted), file(read_encrypted_checksum), file(read_unencrypted_checksum), val(SAMPLE_ID) into read_encrypted

       script:
       read_encrypted = read+".gpg"
       read_encrypted_checksum = read+".gpg.md5"
       read_unencrypted_checksum = read+".md5"
       """
       if [ ! -f ${params.batchDir}/${params.batchName}/${read_encrypted} ]; then
              java -jar ${params.EGAencryptor} -i ${read}
              mv output-files/* .
       else
              ln -s ${params.batchDir}/${params.batchName}/${read_encrypted} ${read_encrypted}
              ln -s ${params.batchDir}/${params.batchName}/${read_encrypted_checksum} ${read_encrypted_checksum}
              ln -s ${params.batchDir}/${params.batchName}/${read_unencrypted_checksum} ${read_unencrypted_checksum}
       fi
       """
}

// Upload files to EGA
process ega_upload {
       cpus 2
       
       input:
       set file(read_encrypted), file(read_encrypted_checksum), file(read_unencrypted_checksum), val(SAMPLE_ID) from read_encrypted

       script:
       path = params.batchName+"/"+SAMPLE_ID
       """
       lftp ftp.ega.ebi.ac.uk -e "cd ${path} || mkdir -p ${path} || cd ${path}; put ${read_encrypted_checksum}; put ${read_unencrypted_checksum}; put ${read_encrypted}; exit"
       touch ${params.batchDir}/${params.batchName}/list_of_ega_uploaded/${read_encrypted}
       """
}

/*
*/