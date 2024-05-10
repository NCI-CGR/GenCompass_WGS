### Visualization of samtools coverage output for multiple samples

## Metric definitions for samtools coverage: http://www.htslib.org/doc/samtools-coverage.html

## Usage
# Rscript plot_samtools_coverage_multi_sample.R "globpath"

## Usage with example data:
# Rscript plot_samtools_coverage_multi_sample.R "/data/COVID_WGS/test_qc_tools/interspecies_contamination/cow_pig_mapping_test/test_mapping_20221201/human_samtools_coverage/*/*.samtoolsCoverage.txt"


### Set up packages and data

library(ggplot2)
library(ggpubr)
#library(dplyr)

theme_set(theme_bw())

# Get filepaths
args = commandArgs(trailingOnly=TRUE)
globpath <- args[1]
filepaths <- Sys.glob(globpath)

# Read in samtools coverage files, create new column with sample name, and combine into a single dataframe
cov <- do.call(rbind, lapply(filepaths, function(filepath){
  sampleName <- gsub(".samtoolsCoverage.txt", "", basename(filepath))
  f <- read.table(filepath, head=T, stringsAsFactors = F, comment.char="")
  f$Sample <- sampleName
  return(f)
  }))

# Fix first column name
names(cov)[1] <- "rname"


### Plot several metrics for main chromosomes

## Subset to main chromosomes and fix order by setting factor levels
chromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))
cov_chrom <- subset(cov, rname %in% chromosomes)
cov_chrom$rname <- factor(cov_chrom$rname, levels=chromosomes)

## For smaller sample sets (n<=10):

if (length(unique(cov$Sample))<=10) {

  ## Plot number of reads per chromosome
  p1 <- ggplot(cov_chrom, aes(x=rname, y=numreads/1e6, color=Sample)) +
    xlab("Chromosome") +
    ylab("Number of Reads (Millions)") +
    scale_x_discrete(labels=c(1:22, "X", "Y", "M")) +
    geom_jitter() + 
    theme(legend.text = element_text(size=8),
          legend.key.size = unit(0.5, 'cm'))
  
  ## Plot % of bases covered per chromosome
  p2 <- ggplot(cov_chrom, aes(x=rname, y=coverage, color=Sample)) +
    xlab("Chromosome") +
    ylab("% of Bases Covered") +
    scale_x_discrete(labels=c(1:22, "X", "Y", "M")) +
    geom_jitter() + 
    theme(legend.text = element_text(size=8),
          legend.key.size = unit(0.5, 'cm'))
  
  ## Plot mean depth per chromosome
  # Plot is truncated at 100x, since mitochondrial depth tends to be very high and make the plot unreadable.
  p3 <- ggplot(cov_chrom, aes(x=rname, y=meandepth, color=Sample)) +
    xlab("Chromosome") +
    ylab("Mean Depth") +
    scale_x_discrete(labels=c(1:22, "X", "Y", "M")) +
    geom_jitter() + 
    theme(legend.text = element_text(size=8),
          legend.key.size = unit(0.5, 'cm')) +
    ylim(c(0,100))
  
  ## Plot mean mapping quality per chromosome
  p4 <- ggplot(cov_chrom, aes(x=rname, y=meanmapq, color=Sample)) +
    xlab("Chromosome") +
    ylab("Mean Mapping Quality") +
    scale_x_discrete(labels=c(1:22, "X", "Y", "M")) +
    geom_jitter() + 
    theme(legend.text = element_text(size=8),
          legend.key.size = unit(0.5, 'cm'))
  
} else {

## For larger sample sets (n>10):

  ## Plot number of reads per chromosome
  p1 <- ggplot(cov_chrom, aes(x=rname, y=numreads/1e6)) +
    xlab("Chromosome") +
    ylab("Number of Reads (Millions)") +
    scale_x_discrete(labels=c(1:22, "X", "Y", "M")) +
    geom_jitter(size=0.5, alpha=0.5) 

  ## Plot % of bases covered per chromosome
  p2 <- ggplot(cov_chrom, aes(x=rname, y=coverage)) +
    xlab("Chromosome") +
    ylab("% of Bases Covered") +
    scale_x_discrete(labels=c(1:22, "X", "Y", "M")) +
    geom_jitter(size=0.5, alpha=0.5)
  
  ## Plot mean depth per chromosome
  # Plot is truncated at 100x, since mitochondrial depth tends to be very high and make the plot unreadable.
  p3 <- ggplot(cov_chrom, aes(x=rname, y=meandepth)) +
    xlab("Chromosome") +
    ylab("Mean Depth") +
    scale_x_discrete(labels=c(1:22, "X", "Y", "M")) +
    geom_jitter(size=0.5, alpha=0.5) +
    ylim(c(0,100))
  
  ## Plot mean mapping quality per chromosome
  p4 <- ggplot(cov_chrom, aes(x=rname, y=meanmapq)) +
    xlab("Chromosome") +
    ylab("Mean Mapping Quality") +
    scale_x_discrete(labels=c(1:22, "X", "Y", "M")) +
    geom_jitter(size=0.5, alpha=0.5)

}


### Save plots to PDF
plots <- ggarrange(p1, p2, p3, p4, nrow=4)
ggsave("samtools_coverage_plots_by_chromosome.pdf", plots, width=8, height=12)


### Print list of sequences with no reads detected

## Print individual sequence and sample names
# cov[cov$numreads==0, c("Sample", "rname")]  

## Print table (sequence name and number of samples that had 0 reads detected)
tab <- as.data.frame(table(cov[cov$numreads==0, "rname"]))
names(tab) <- c("Sequence_Name", "Num_Samples_0_Reads")
write.table(tab, "sequences_with_zero_reads.txt", quote=F, row.names=F, sep="\t")
