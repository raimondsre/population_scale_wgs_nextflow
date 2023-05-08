args <- commandArgs(trailingOnly = TRUE)
lapply(c("data.table","dplyr","ggplot2","tidyr"), require, character.only = TRUE)
          
input_file <- NULL
for (i in seq_along(args)) {
    if (args[i] == "--input") {input_file <- args[i + 1]}
    else if (args[i] == "--interval") {interval <- args[i + 1]}
    else if (args[i] == "--original_file_name") {original_file_name <- args[i + 1]}
}
# Read input file using fread
if (!is.null(input_file)) {
    data <- fread(input_file, h=TRUE)
} else {stop('input file not specified')}
if (is.null(interval)) { interval <- "all" } 

# count variants in specific INFO groups
data <- data %>% 
mutate(AF_GROUP = ifelse(AF >= 0.05 & AF <= 0.95, "Common_variants", ifelse(AF < 0.005 | AF > 0.995, "Rare_variants", "Low_frequency_variants"))) %>%
mutate(snv = ifelse(nchar(REF) == 1 & nchar(ALT) == 1,TRUE,FALSE)) %>%
mutate(INFO_GROUP = ifelse(INFO >= 0.8, "high_confidence", ifelse(INFO > 0.4 & INFO < 0.8, "low_confidence", "Below0_4_confidence"))) %>%
group_by(AF_GROUP, snv, INFO_GROUP) %>%
summarise(count = n())


# Add interval and sample
total <- data %>% mutate(interval = interval, original_file_name = original_file_name)
fwrite(total, sprintf("%s.counted.txt", original_file_name),col.names=FALSE, sep="\t")