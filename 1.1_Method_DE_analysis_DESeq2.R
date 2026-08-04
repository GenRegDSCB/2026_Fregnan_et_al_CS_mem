## Package install -----
library("tidyverse")
library("DESeq2")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

sdata <- read.delim("0_Input_metadata.tsv", row.names=1)

rnaseq <- read.delim("0_Input_raw_counts.tsv", 
                     row.names=1, check.names = F)

fpkm <- read.delim("0_Input_FPKM.tsv", 
                   row.names=1, check.names = F)

sdata$Class <- ifelse(sdata$Class == "PENIS_CS_FLAT", "PENIS_CS", sdata$Class)
sdata$Class <- ifelse(sdata$Class == "PENIS_CS_GRAT", "PENIS_CS", sdata$Class)
sdata$Class <- ifelse(sdata$Class == "PENIS_CS_SCA", "PENIS_CS", sdata$Class)

## DE analysis preprocessing -----
countdata <- rnaseq[, row.names(sdata)]
fpkm <- fpkm[, row.names(sdata)]

countdata <- round(countdata)

dds_complete <- DESeqDataSetFromMatrix(countData = countdata, 
                                       colData = sdata, 
                                       design = ~Class)

dds_complete <- DESeq(dds_complete)

### Selection of the classes -----
treatment <- "CN_CS_MEM_SCA"
control <- "CN_CS_MEM_FLAT"

## DE analysis testing -----
res_sub <- results(dds_complete,
                   contrast=c("Class", treatment, control),
                   cooksCutoff=FALSE)

res_sub <- data.frame(res_sub) %>% rownames_to_column('ID')

write.table(res_sub, 
            paste0("2_DE_analysis_results_", treatment, "_", control, ".tsv"), 
            sep="\t", 
            quote=F, 
            row.names=F)

### Median levels ----
norm_counts <- counts(dds_complete, normalized=T)
write.table(norm_counts, "0_Input_normalized_counts.tsv", sep="\t", quote=F)

dat <- data.frame(t(fpkm), Class = sdata$Class, check.names = F)
dat <- gather(dat, Gene, Level, -Class)
dat <- group_by(dat, Gene, Class) %>% dplyr::summarize(med=median(Level))
out <- pivot_wider(dat, names_from = Class, values_from = med)

write.table(out, "2.2_Table_Median_counts.tsv", sep="\t", quote=F, row.names=F)