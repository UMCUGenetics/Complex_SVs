# This script is used to collect and combine the bamMetrics files for merged BAM files. 
# These bamMetrics files are used for quality control

args <- commandArgs(trailingOnly=TRUE)
filelist <- args[1]
output_folder <- args[2]

dfilelist <- read.delim(filelist, header = F, stringsAsFactors = F)

for(i in 1:nrow(dfilelist)){
  sample <- dfilelist$V1[i]
  path1 <- dfilelist$V2[i]
  path2 <- dfilelist$V3[i]
  
  readaligmnet1 <- read.delim(paste(path1, "/bamMetrics/", sample, "_sort/", sample, "_sort_MultipleMetrics.txt.alignment_summary_metrics", sep = ""), skip = 6 )
  readaligmnet2 <- read.delim(paste(path2, "/bamMetrics/", sample, "_sort/", sample, "_sort_MultipleMetrics.txt.alignment_summary_metrics", sep = ""), skip = 6 )
  
  data <- NULL
  for(i in 1:nrow(readaligmnet1)){
    TOTAL_READS <- sum(readaligmnet1$TOTAL_READS[i], readaligmnet2$TOTAL_READS[i])
    PF_HQ_ALIGNED_READS <- sum(readaligmnet1$PF_HQ_ALIGNED_READS[i], readaligmnet2$PF_HQ_ALIGNED_READS[i])
    PF_READS_ALIGNED <- sum(readaligmnet1$PF_READS_ALIGNED[i], readaligmnet2$PF_READS_ALIGNED[i])
    
    data <- rbind(data, cbind(TOTAL_READS, PF_HQ_ALIGNED_READS, PF_READS_ALIGNED))
  }
  
  print(sample)
  print(data)
  readaligmnet3 <- setNames(data.frame(matrix(ncol = length(colnames(readaligmnet1)), nrow = 3)), colnames(readaligmnet1))
  readaligmnet3[,colnames(data)] <- data
  readaligmnet3[is.na(readaligmnet3)] <- 0
  
  print(readaligmnet3)
  
  dir.create(paste(output_folder, "/bamMetrics/", sep = ""), showWarnings = F)
  outputfile <- paste(output_folder, "/bamMetrics/", sample, "_sort/", sample, "_sort_MultipleMetrics.txt.alignment_summary_metrics", sep = "")
  dir.create(paste(output_folder, "/bamMetrics/", sample, "_sort/", sep = ""), showWarnings = F)
  header <- data.frame(c(7:2))
  write.table(header, outputfile, quote = F, row.names = F, col.names = F)
  write.table(readaligmnet3, outputfile, quote = F, row.names = F, col.names = T, sep = "\t", append = T)

}



RNAMS_all <- NULL
for(i in 1:nrow(dfilelist)){
  RNAMS3 <- data.frame(matrix(nrow = 1))
  sample <- dfilelist$V1[i]
  path1 <- dfilelist$V2[i]
  path2 <- dfilelist$V3[i]
  
  RNAMS1 <- read.delim(paste(path1, "/bamMetrics/", "RNAMetrics_summary.txt", sep = ""), header = T)
  RNAMS1 <- RNAMS1[grep(sample, RNAMS1$sample),][,1:11]
  
  RNAMS2 <- read.delim(paste(path2, "/bamMetrics/", "RNAMetrics_summary.txt", sep = ""), header = T)
  RNAMS2 <- RNAMS2[grep(sample, RNAMS2$sample),][,1:11]
  
  
  RNAMS3$PF_BASES <- sum(RNAMS1$PF_BASES, RNAMS2$PF_BASES)
  RNAMS3$PF_ALIGNED_BASES <- sum(RNAMS1$PF_ALIGNED_BASES, RNAMS2$PF_ALIGNED_BASES)  
  RNAMS3$RIBOSOMAL_BASES <- sum(as.numeric(RNAMS1$RIBOSOMAL_BASES), as.numeric(RNAMS2$RIBOSOMAL_BASES))
  RNAMS3$CODING_BASES <- sum(RNAMS1$CODING_BASES, RNAMS2$CODING_BASES)
  RNAMS3$UTR_BASES <- sum(RNAMS1$UTR_BASES, RNAMS2$UTR_BASES)
  RNAMS3$INTRONIC_BASES <- sum(RNAMS1$INTRONIC_BASES, RNAMS2$INTRONIC_BASES)
  RNAMS3$INTERGENIC_BASES <- sum(RNAMS1$INTERGENIC_BASES, RNAMS2$INTERGENIC_BASES)
  RNAMS3$IGNORED_READS <- sum(RNAMS1$IGNORED_READS, RNAMS2$IGNORED_READS)
  RNAMS3$CORRECT_STRAND_READS <- sum(RNAMS1$CORRECT_STRAND_READS, RNAMS2$CORRECT_STRAND_READS)
  RNAMS3$INCORRECT_STRAND_READS <- sum(RNAMS1$INCORRECT_STRAND_READS, RNAMS2$INCORRECT_STRAND_READS)
  RNAMS3$matrix.nrow...1. <- paste(sample, "_sort", sep = "")
  colnames(RNAMS3) <- colnames(RNAMS2)
  
  RNAMS3$PCT_RIBOSOMAL_BASES <- RNAMS3$RIBOSOMAL_BASES/RNAMS3$PF_ALIGNED_BASES
  RNAMS3$PCT_CODING_BASES <- RNAMS3$CODING_BASES/RNAMS3$PF_ALIGNED_BASES
  RNAMS3$PCT_UTR_BASES <- RNAMS3$UTR_BASES/RNAMS3$PF_ALIGNED_BASES
  RNAMS3$PCT_INTRONIC_BASES <- RNAMS3$INTRONIC_BASES/RNAMS3$PF_ALIGNED_BASES
  RNAMS3$PCT_INTERGENIC_BASES <- RNAMS3$INTERGENIC_BASES/RNAMS3$PF_ALIGNED_BASES
  RNAMS3$PCT_MRNA_BASES <- NA
  RNAMS3$PCT_USABLE_BASES <- NA
  RNAMS3$PCT_CORRECT_STRAND_READS <- NA
  RNAMS3$MEDIAN_CV_COVERAGE  <- NA 
  RNAMS3$MEDIAN_5PRIME_BIAS <- NA
  RNAMS3$MEDIAN_3PRIME_BIAS  <- NA
  RNAMS3$MEDIAN_5PRIME_TO_3PRIME_BIAS  <- NA 
  RNAMS3$SAMPLE <- NA 
  RNAMS3$LIBRARY <- NA
  RNAMS3$READ_GROUP <- NA
  RNAMS_all <- rbind(RNAMS_all, RNAMS3)
}

outputfile <- paste(output_folder, "/bamMetrics/RNAMetrics_summary.txt", sep = "")
write.table(RNAMS_all, outputfile, quote = F, row.names = F, col.names = T, sep = "\t")




