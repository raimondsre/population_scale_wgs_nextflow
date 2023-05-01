args <- commandArgs(trailingOnly = TRUE)
lapply(c("data.table","dplyr","ggplot2","tidyr"), require, character.only = TRUE)
          
input_file <- NULL
for (i in seq_along(args)) {
    if (args[i] == "--input") {input_file <- args[i + 1]}
    else if (args[i] == "--interval") {interval <- args[i + 1]}
    else if (args[i] == "--sample") {sample <- args[i + 1]}
    else if (args[i] == "--original_file_name") {original_file_name <- args[i + 1]}
}
# Read input file using fread
if (!is.null(input_file)) {
    data <- fread(input_file, h=FALSE)
} else {stop('VEP annotated input file not specified')}
if (is.null(interval)) { interval <- "all" } 
if (is.null(sample)) { sample <- "all" }  
# If annotated file empty, write empty file and stop
if (length(data) == 0) {
    file.create(sprintf("%s.%s.vep.counted", interval, sample))
    print('VEP annotated input file empty')
    q("yes")
}

# Separate into individual features
vep <- data %>% 
    setNames(c("snpID","var_type","clinical","position","rsID","gnomAD_freq","associated_phenotype")) %>%
    separate_rows(., clinical, sep = "&") %>%
    separate_rows(., position, sep = "&") %>%
    separate_rows(., clinical, sep = "/") %>%
    separate_rows(., position, sep = "/") %>%
    mutate(novel = ifelse(rsID == ".",TRUE,FALSE))
# Count occurance of each unique feature
total <- rbind(
    # Existing
    vep %>% filter(!duplicated(paste(snpID,var_type))) %>% count(var_type, name = "Count") %>% setNames(c("var","count")) %>% mutate(novel_vep=1),
    vep %>% filter(!duplicated(paste(snpID,clinical))) %>% count(clinical, name = "Count") %>% setNames(c("var","count")) %>% mutate(novel_vep=1),
    vep %>% filter(!duplicated(paste(snpID,position))) %>% count(position, name = "Count") %>% setNames(c("var","count")) %>% mutate(novel_vep=1),
    # Novel
    vep[vep$novel,] %>% filter(!duplicated(paste(snpID,var_type))) %>% count(var_type, name = "Count") %>% setNames(c("var","count")) %>% mutate(novel_vep=0),
    vep[vep$novel,] %>% filter(!duplicated(paste(snpID,clinical))) %>% count(clinical, name = "Count") %>% setNames(c("var","count")) %>% mutate(novel_vep=0),
    vep[vep$novel,] %>% filter(!duplicated(paste(snpID,position))) %>% count(position, name = "Count") %>% setNames(c("var","count")) %>% mutate(novel_vep=0)
)
# Add interval and sample
total <- total %>% mutate(interval = interval, sample = sample, original_file_name = original_file_name)
fwrite(total, sprintf("%s.%s.vep.counted", interval, sample),col.names=FALSE, sep="\t")