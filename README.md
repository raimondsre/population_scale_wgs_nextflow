# population_scale_wgs_nextflow
Nextflow pipelines for population scale WGS analysis

Complete workflow replicating main results of Latvian Genome Reference paper will be produced running command
nextflow run raimondsre/population_scale_wgs_nextflow -r main --input reference_panel.vcf.gz

Subworkflows can be run separately, e.g. 
nextflow run raimondsre/population_scale_wgs_nextflow/templates/imputation_panel_comparison.nf -r main --input reference_panel.vcf.gz