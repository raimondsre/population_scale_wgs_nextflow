#!/usr/bin/env nextflow

params.sampleLocation = './for.trimming'
params.batchDir = projectDir
params.batchName = 'lv_reference_20220722_502samples'
params.EGAencryptor = '/home/raimondsre/programms/EGA-Cryptor-2.0.0/ega-cryptor-2.0.0.jar'

def remPath(String fileName) {return fileName.replaceAll(/.*\//,'')}

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
       .map {tuple(it, (params.batchDir+"/"+params.batchName+"/"+remPath(it[1])+".gpg").exists)}
       .into { test }
test.subscribe {println it}
/*
       .filter { (params.batchDir+"/"+params.batchName+"/"+remPath(it[1])+".gpg").isEmpty() }
       .into { to_encrypt; intervals2 }
to_encrypt.subscribe { println it}

/*
// Encrypt each file
process ega_encrypt {
       publishDir params.batchDir+"/"+params.batchName, mode: 'copy'

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
       java -jar ${params.EGAencryptor} -i ${read}
       mv output-files/* .
       """
}

// Upload files to EGA

/*
*/