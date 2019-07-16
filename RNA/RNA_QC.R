## This QC is not essential for the differential expression analysis.
# The number of reads per sample and the mapping (for example how many reads map to exons) are checked

library(ggplot2)
library(reshape2)
library(gridExtra)

options(scipen = 999)

args <- commandArgs(trailingOnly=TRUE)
filelist_path <- args[1]
output_folder <- args[2]
output_name <- args[3]

output_file <- paste(output_name, "_RNA_QC_Report", sep = "")


out <- paste(output_folder, "QC/", sep = "")
idxstats_folder <- paste(out, "Idxstats/", sep = "")

ifelse(!dir.exists(file.path(out)), dir.create(file.path(out), showWarnings = TRUE), FALSE)

ifelse(!dir.exists(file.path(idxstats_folder)), dir.create(file.path(idxstats_folder), showWarnings = TRUE), FALSE)


print("### Starting RNA-seq Quality Control ###")
print(paste("# Output folder = ", out, sep = ""))
# The filelist is used to determine the locations of the BAM files.
# The filelist should contain the following columns: 
#1 Subject (this is the patient/subject ID), 
#2 Sample (this is the sample ID, some patients have multiple samples which will be used as replicates in differential expression), 
#3 Folder (Folder ), 
#4 Run (This is the ID for the sequencing run, for quality comparison between runs), 
#5 Sex (Gender of the subject, for checking X/Y reads)
filelist <- read.delim(filelist_path, header = T, stringsAsFactors = F)

#print(paste("Read filelist. Number of samples = ", nrow(filelist), ". Head filelist: ", sep =""))
#print(head(filelist))

plot_qc_data <- function(filelist, relative = FALSE, ylab = "Reads (millions)", title = "Aligned PF reads", write = TRUE, png = TRUE){
  alignment_overview <- NULL
  
  for(i in 1:nrow(filelist)){
  	sample <- filelist$Sample[i]
    folder <- filelist$Folder[i]
    
    alignment_summary_file <- paste(folder, "/bamMetrics/", sample, "_sort/", sample, "_sort_MultipleMetrics.txt.alignment_summary_metrics", sep = "")
    
    alignment_summary <- read.delim(alignment_summary_file, skip = 6)
    alignment_summary$sample <- sample
    alignment_summary$run <- filelist$Run[i]
    alignment_overview <- rbind(alignment_overview, alignment_summary[3,])
    
  }  
  
  if(write == TRUE){
    write.table(alignment_overview, file = paste(output_folder, format(Sys.Date(), "%Y%m%d"), "_RNA_alignment_QC.txt", sep = ""), sep = "\t")
  }
  
  qc_plot_data <- data.frame(sample = alignment_overview$sample,
                             run = alignment_overview$run,
                             aligned_reads_Q20 = alignment_overview$PF_HQ_ALIGNED_READS,
                             aligned_reads = alignment_overview$PF_READS_ALIGNED - alignment_overview$PF_HQ_ALIGNED_READS,
                             unaligned_reads = alignment_overview$TOTAL_READS - alignment_overview$PF_READS_ALIGNED)
  
  average <- data.frame(sample = "Average", 
                        run = "Average",
                        aligned_reads_Q20 = as.numeric(mean(qc_plot_data[,"aligned_reads_Q20"])),
                        aligned_reads =as.numeric(mean(qc_plot_data[,"aligned_reads"])),
                        unaligned_reads = as.numeric(mean(qc_plot_data[,"unaligned_reads"])))
  
  
  qc_plot_data <- rbind(qc_plot_data, average)

  if(relative == TRUE){
    division <- 1
    for(row in 1:nrow(qc_plot_data)){
      qc_plot_data[row,-1] <- qc_plot_data[row,-1]/sum(qc_plot_data[row,-1]) * 100
      ylab <- "% of Reads"
      title <- "Aligned PF reads (relative)"
    }
  } else {
    division <- 1e6
  }
  
  qc_plot_data_m <- melt(qc_plot_data)     
  qc_plot_data_m$variable <- factor(qc_plot_data_m$variable, levels = c("unaligned_reads","aligned_reads","aligned_reads_Q20" ))
  qc_plot_data_m$color <- factor(qc_plot_data_m$variable, levels = c("unaligned_reads","aligned_reads","aligned_reads_Q20" ))
  #qc_plot_data_m$sample <- factor(qc_plot_data_m$sample, levels = qc_plot_data_m$sample[order(qc_plot_data_m$run)])
  levels(qc_plot_data_m$color) <- c("#F8766D","#619CFF","#00BA37") 
  
  qc_plot <- ggplot(qc_plot_data_m, aes(x = sample, y = value/division, fill = variable)) + 
    geom_bar(width = 0.7,stat = "identity") + 
    scale_fill_manual(values = levels(qc_plot_data_m$color)) + 
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Sample", y = ylab) +
    ggtitle(title) +
    geom_abline(intercept = qc_plot_data[qc_plot_data$sample == "Average","aligned_reads_Q20"]/division, slope = 0) + 
    geom_abline(intercept = sum(qc_plot_data[qc_plot_data$sample == "Average",c("aligned_reads_Q20","aligned_reads")])/division, slope = 0)+
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_grid(~run, scales = "free", space = "free") +
    theme(panel.grid.major.y = element_line(colour = "gray75"), panel.grid.minor.y = element_line(colour = "gray90"))
  qc_plot
}


plot_annotation_data <- function(filelist, ylab = "% of bases", title = "RNA Metrics", relative = TRUE, write = TRUE, png = TRUE){
  print("# Plotting annotation data")
  RNAMetrics_overview <- NULL
  for(folder in unique(filelist$Folder)){
    #print(folder)
    run <- unique(filelist[filelist$Folder == folder, "Run"])
    RNAMetrics <- read.delim(paste(folder, "/bamMetrics/RNAMetrics_summary.txt", sep = ""))
    RNAMetrics$run <- run
    RNAMetrics_overview <- rbind(RNAMetrics_overview,RNAMetrics)
    
    
    RNAMetrics_overview$sample <- gsub(pattern = "_sort", replacement = "", RNAMetrics_overview$sample)
    RNAMetrics_overview <- RNAMetrics_overview[RNAMetrics_overview$sample %in% filelist$Sample,]
  }
  
  if(write == TRUE){
    write.table(RNAMetrics_overview, file = paste(out,output_name, "_RNAMetrics_QC.txt", sep = ""), sep = "\t")
  }
    
  RNAMetrics_to_plot <- data.frame(sample = RNAMetrics_overview$sample,
                                   run = RNAMetrics_overview$run,
                                   Ribosomal = RNAMetrics_overview$PCT_RIBOSOMAL_BASES,
                                   Coding = RNAMetrics_overview$PCT_CODING_BASES,
                                   UTR = RNAMetrics_overview$PCT_UTR_BASES,
                                   Intronic = RNAMetrics_overview$PCT_INTRONIC_BASES,
                                   Intergenic = RNAMetrics_overview$PCT_INTERGENIC_BASES)
  
  average <- data.frame(sample = "Average", 
                        run = "Average",
                        Ribosomal = mean(RNAMetrics_to_plot[,"Ribosomal"]),
                        Coding = mean(RNAMetrics_to_plot[,"Coding"]),
                        UTR = mean(RNAMetrics_to_plot[,"UTR"]),
                        Intronic = mean(RNAMetrics_to_plot[,"Intronic"]),
                        Intergenic = mean(RNAMetrics_to_plot[,"Intergenic"]))
  RNAMetrics_to_plot <- rbind(RNAMetrics_to_plot, average)
  

  if(relative == TRUE){
    RNAMetrics_to_plot_m <- melt(RNAMetrics_to_plot)
    RNAMetrics_to_plot_m$variable <- factor(RNAMetrics_to_plot_m$variable, levels = c("Ribosomal", "Intergenic","Intronic", "UTR","Coding"))
    
    
    RNAMetrics_plot <- ggplot(RNAMetrics_to_plot_m, aes(x = sample, y = value*100, fill = variable)) + 
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = c("Ribosomal" = "#F8766D", "Intergenic" = "#E76BF3","Intronic" = "#A3A500","UTR" ="#619CFF", "Coding"="#00BA37"))  + 
      facet_grid(~run, scales = "free", space = "free")+
      theme_minimal()+
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(x = "Sample", y = ylab)+
      theme(axis.text.x = element_text(angle = 90)) +
      geom_abline(intercept = sum(RNAMetrics_to_plot[RNAMetrics_to_plot$sample == "Average",c("Coding", "UTR")])*100, slope = 0) + 
      theme(panel.grid.major.y = element_line(colour = "gray75"), panel.grid.minor.y = element_line(colour = "gray90"))
    RNAMetrics_plot
  } else if (relative == FALSE){
    
    Usuable_reads <- data.frame(sample = RNAMetrics_overview$sample, Coding = RNAMetrics_overview$CODING_BASES/75, UTR = RNAMetrics_overview$UTR_BASES/75, Run = RNAMetrics_overview$run)
    Usuable_reads_m <- melt(Usuable_reads)
    ggplot(Usuable_reads_m, aes(x = sample, fill = variable, y = value)) + geom_bar(stat = "identity", width = 0.7) +  facet_grid(~Run, scales = "free", space = "free") +
      theme_minimal()+ theme(axis.text.x = element_text(angle = 90)) + labs(y = "Reads (Number)") +
      scale_fill_manual(values = c("UTR" ="#619CFF", "Coding"="#00BA37"))
    
  }

}

run_idxstats <- function(filelist){
  print("# Running Samtools idxstats")
  for(i in 1:nrow(filelist)){
    #print(i)
    sample <- filelist[i,"Sample"]
    folder <- filelist[i,"Folder"]
    
    bam_folder <- paste(folder, "/", sample, "/mapping/", sep = "")
    
    bam_file  <- list.files(bam_folder, pattern = "\\.bam$")
    
    bam <- paste(bam_folder, bam_file, sep = "")
    
    command <- paste("samtools idxstats ", bam, " > ", idxstats_folder, sample, "_idxstats.txt", sep  ="")
    
    print(command)
    system(command)
  }
  
}

gender_check <- function(filelist, ylab = "Reads (millions)"){
  print("# Running Gender check")
  all_stats <- data.frame(row.names = c(1:22, "X", "Y", "MT", "*"))
  for(i in 1:nrow(filelist)){
    #print(i)
    sample <- filelist[i,"Sample"]
    folder <- filelist[i,"Folder"]
    run <- filelist[i,"Run"]
    
    gender <- filelist[i,"Sex"]
    
    idxstats <- read.delim(paste(idxstats_folder,sample, "_idxstats.txt",sep = ""), header = F, row.names = 1)
    
    idxstats["*",2] <- idxstats["*",3]
    all_stats <- cbind(all_stats, idxstats[,2])
    names(all_stats)[length(names(all_stats))] <- as.character(sample)
  }
  
  rownames(all_stats)[rownames(all_stats) == "*"] <- "Unaligned"
  all_stats$Average <- rowMeans(all_stats)
  
  males <- filelist[filelist$Sex == "M", "Sample"]
  females <- filelist[filelist$Sex == "F", "Sample"]
  
  all_stats$Average_males <- rowMeans(all_stats[,names(all_stats) %in% males])
  all_stats$Average_females <- rowMeans(all_stats[,names(all_stats) %in% females])
  all_stats["Total",] <- colSums(all_stats)
  
  ## Gender-check
  all_stats_sex <- all_stats[c("Total", "Y"),] # select data only for X and Y chromosome
  all_stats_sex["Y_rel",] <- all_stats_sex["Y",] / all_stats_sex["Total",] # divide #reads Y-chromosome by # reads X-chromosome
  
  all_stats_sex_m <- melt(all_stats_sex["Y_rel",])
  
  for(i in 1:nrow(all_stats_sex_m)){
    if( as.character(all_stats_sex_m[i,"variable"]) %in% filelist$Sample){
      all_stats_sex_m[i,"sex"] <- as.character(filelist[filelist$Sample == as.character(all_stats_sex_m[i,"variable"]), "Sex"])
      all_stats_sex_m[i,"run"] <- as.character(filelist[filelist$Sample == as.character(all_stats_sex_m[i,"variable"]), "Run"])
      #all_stats_sex_m[i,"prep"] <- as.character(filelist[filelist$Sample == as.character(all_stats_sex_m[i,"variable"]), "Prep"])
    }
  }  
  gender_plot <- ggplot(data = all_stats_sex_m, aes(x = variable, y = value*100, fill = sex)) +
    geom_bar(stat = "identity") + 
    geom_text(data = all_stats_sex_m, aes(x = variable, y = value*100+0.05, label = sex)) +
    geom_abline(intercept = all_stats_sex_m[all_stats_sex_m$variable == "Average_males","value"]*100, slope = 0, col = "darkblue") + 
    geom_abline(intercept = all_stats_sex_m[all_stats_sex_m$variable == "Average_females","value"]*100, slope = 0, col = "purple") + 
    facet_grid(~run,scales = "free", space = "free") +
    ggtitle("Gender Check") + labs(y = "% of reads on chrY") + theme_minimal() + theme(axis.text.x = element_text(angle = 90)) + 
    theme(panel.grid.major = element_line(colour = "gray75"), panel.grid.minor = element_line(colour = "gray90")) +
    theme(plot.title = element_text(hjust = 0.5))
  gender_plot
  
}






# ######## COUNT QC

run_idxstats(filelist = filelist)

pdf(file = paste(out,output_file, ".pdf", sep= ""), pointsize = 8, width = 11, height = 8)

plot_qc_data(filelist = filelist, relative = FALSE)
plot_annotation_data(filelist = filelist, relative = FALSE)
plot_annotation_data(filelist = filelist, relative = TRUE)

gender_check(filelist = filelist)


dev.off()

png(paste(out, output_name, "_RNAMetrics_QC.png", sep= ""), width = 1280, height = 720)
plot_annotation_data(filelist = filelist)
dev.off()

png(paste(out, output_name, "_RNA_alignment_QC.png", sep= ""), width = 1280, height = 720)
plot_qc_data(filelist = filelist, relative = FALSE)
dev.off()


print("### Finished RNA-seq Quality Control ###")