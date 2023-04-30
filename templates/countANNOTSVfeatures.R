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
    data <- fread(input_file, h=TRUE)
} else {stop('VEP annotated input file not specified')}
if (is.null(interval)) { interval <- "all" } 
if (is.null(sample)) { sample <- "all" }  


annotsv <- data %>% 
    filter(!duplicated(AnnotSV_ID)) %>%
    select(AnnotSV_ID,SV_type,GnomAD_pLI,B_loss_source,B_gain_source) %>% 
    mutate(novel_annotsv = ifelse(is.na(GnomAD_pLI) & 
        B_loss_source=="" & 
        B_gain_source=="",
        "TRUE", "FALSE")) %>%
    group_by(novel_annotsv,SV_type) %>%
    count(SV_type, name = "Count") %>% 
    select(var=SV_type,count=Count,novel_annotsv) %>%
    mutate(interval = interval,
        sample = sample,
        original_file_name = original_file_name) %>%
    fwrite(sprintf("%s.%s.annotsv.counted", interval, sample), col.names=FALSE, sep="\t")