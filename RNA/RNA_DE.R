# This scripts calculates the differential expression of all subjects in the filelist.
# It also does some QC including clustering of samples. 
# This script can be started using the RNA_DE.sh script, else the arguments will have to be added manually

library(DESeq2)
library(ggplot2)
library(gplots)
library(reshape2)
library(RColorBrewer)
library(genefilter)
library(edgeR)
library(biomaRt)

args <- commandArgs(trailingOnly=TRUE)
filelist_path <- args[1]
output_folder <- args[2]
output_name <- args[3]
qc <- args[4]
plots <- args[4]
truncated_genes <- args[5]
GTFFile_TruncatedGenes <- args[6]
Data_TruncatedGenes <- args[7]
Exon_sizes <- args[8]

# The filelist is used to determine the locations of the BAM files.
# The filelist should contain the following columns: 
#1 Subject (this is the patient/subject ID), 
#2 Sample (this is the sample ID, some patients have multiple samples which will be used as replicates in differential expression), 
#3 Folder (Folder ), 
#4 Run (This is the ID for the sequencing run, for quality comparison between runs), 
#5 Sex (Gender of the subject, for checking X/Y reads)
#6 FolderTruncatedGenes (Folder to HTseqCount files of the truncated genes)
filelist <- read.delim(filelist_path, header = TRUE, stringsAsFactors = FALSE)

options(scipen = 999)

DE_output <- paste(output_folder, "DE/", sep = "")

ifelse(!dir.exists(file.path(DE_output)), dir.create(file.path(DE_output), showWarnings = FALSE), FALSE)

# First the raw RNA counts (HTSeq output) and RPKMs for each sample are collected and merged to the rna_counts dataframe.

if(truncated_genes == FALSE){
  print("Merging raw count tables")
  rna_counts <- data.frame()
  rna_rpkms <- data.frame()
  
  # collect the raw_counts and rpkm files from each sequencing run (which are located in different folders)
  for(input_folder in unique(filelist$Folder)){
    print(input_folder)

    count_file <- paste(input_folder, "/read_counts/", list.files(path = paste(input_folder,"/read_counts", sep = ""), pattern = "readCounts_raw.txt"), sep = "")
    counts <- read.delim(file = count_file, header = T, stringsAsFactors = F, check.names = F)
  
    counts <- counts[,names(counts) %in% filelist$Sample[filelist$Folder == input_folder] | names(counts) == "gene"] # only select the columns (=samples) that are included in the genelist
    # (there may be samples in the counts file that should not be included in DE analysis)

    # Add the counts from this count file to the counts of other folder/runs
    if(nrow(rna_counts) == 0){
        rna_counts <- counts
      } else { 
        rna_counts <- merge(rna_counts, counts, by = "gene")
      }
    
    rpkms_file <- paste(input_folder, "/read_counts/", list.files(path = paste(input_folder,"/read_counts", sep = ""), pattern = "RPKM.txt"), sep = "")
    rpkms <- read.delim(rpkms_file, header = T, stringsAsFactors = F, check.names = F)
    names(rpkms) <- gsub("[.]", "-", names(rpkms)) #The samples names in the RPKM file contain "." instead of "-" 
    names(rpkms)[1] <- "gene" # first column in RPKM file contains ensembl ids, but it does not have a name
    
    # Some sample names (columns) in the RPKM files start with a "X". Remove this X if present:
    names(rpkms)[startsWith(names(rpkms), "X")] <- substring(names(rpkms)[startsWith(names(rpkms), "X")], 2)
  
    rpkms <- rpkms[,names(rpkms) %in% filelist$Sample[filelist$Folder == input_folder] | names(rpkms) == "gene"] # only select the columns (=samples) that are included in the genelist 
    # (there may be samples in the counts file that should not be included in DE analysis)
  
    if(nrow(rna_rpkms) < 1){
      rna_rpkms <- rpkms
    } else { 
      rna_rpkms <- merge(rna_rpkms, rpkms, by = "gene")
    }
  }

  # Set the gene names as row.names and remove the "gene" column
  row.names(rna_counts) <- rna_counts$gene
  rna_counts <- rna_counts[,-which(names(rna_counts) == "gene")]
  
  row.names(rna_rpkms) <- rna_rpkms$gene
  rna_rpkms <- rna_rpkms[,-which(names(rna_rpkms) == "gene")]

  print(paste("Writing merged raw RNA counts table to: ", output_folder, output_name, "_RNA_raw_counts.txt", sep = ""))
  write.table(rna_counts, file = paste(output_folder, output_name, "_RNA_raw_counts.txt", sep = ""),quote = F, sep = "\t", row.names = T)
  write.table(rna_rpkms, file = paste(output_folder, output_name, "_RNA_RPKMs.txt", sep = ""),quote = F, sep = "\t", row.names = T)
  print(paste("Merged raw RNA count tables of ", ncol(rna_counts), " samples", sep =""))
}

if(truncated_genes == TRUE){
  print("Merging raw count tables with truncated gene fragements included")
  
  
  ## Determine which truncated fragements are not counted > intronic fragments
  gtf_trunc <- read.delim(file = GTFFile_TruncatedGenes, header = F, stringsAsFactors = F)
  present_ensembl_trunc <- NULL
  for(i in 1:nrow(gtf_trunc)){
    ensembl_id <- gsub(";.*", "", gtf_trunc$V9[i])
    ensembl_id <- gsub("gene_id ", "", ensembl_id)
    present_ensembl_trunc <- rbind(present_ensembl_trunc, ensembl_id)
  }
  present_ensembl_trunc <- as.data.frame(unique(present_ensembl_trunc))
  d_TruncatedGenes <- read.delim(Data_TruncatedGenes, header = T, stringsAsFactors = F)
  d_frag <- paste(paste(d_TruncatedGenes$ensembl_gene_id, "_", sep = ""), d_TruncatedGenes$Fragment, sep = "")
  absent_frag <- as.data.frame(d_frag[!d_frag %in% present_ensembl_trunc$V1])
  absent_frag$V2 <- 0
  
  rna_counts <- data.frame()
  for(i in 1:nrow(filelist)){
    input_folder <- filelist$Folder[i]
    input_folderTruncatedGenes <- filelist$FolderTruncatedGenes[i]
    Sample <- filelist$Sample[i]
    
    # HTSeq counts complete genes
    count_file <- paste(input_folder, "/", Sample, "/read_counts/", Sample, "_htseq_counts.txt", sep = "")  ###
    counts <- read.delim(file = count_file, header = F, stringsAsFactors = F, check.names = F)
    colnames(counts) <- (c("gene", Sample))
    
    # HTSeq counts truncated genes
    count_file_TG <- paste(input_folderTruncatedGenes, "/", Sample, "/", output_name, "_", Sample, "_HTSeq_Truncated_count.txt", sep = "")   
    if (file.exists(count_file_TG)){
      counts_TG <- read.delim(file = count_file_TG, header = F, stringsAsFactors = F, check.names = F)
      print(count_file_TG)
    }else{
      print(paste("No truncated gene file for", Sample, sep = " "))
      next
    }
      colnames(counts_TG) <- (c("gene", Sample))
    
    # Remove not counted reads
    counts_TG <- counts_TG[grep("ENSG", counts_TG$gene),]
    
    # Add absent fragments from GTFFile 
    colnames(absent_frag) <- c("gene", Sample)
    counts_TG <- rbind(counts_TG, absent_frag)
    
    # Combine counts complete and truncated genes
    counts <- rbind(counts_TG, counts)
    
    if(nrow(rna_counts) < 1){
      rna_counts <- counts
    } else { 
      rna_counts <- merge(rna_counts, counts, by = "gene")
    }
  }
  
  # Remove intact genes of the truncated variants
  intact_genes <- strsplit(counts_TG$gene, "_", fixed = T)
  intact_genes <- data.frame(do.call(rbind, intact_genes))
  rna_counts <- rna_counts[!(rna_counts$gene %in% intact_genes$X1),]

  # Set the gene names as row.names
  row.names(rna_counts) <- rna_counts$gene
  rna_counts <- rna_counts[,-1]
  
  print(paste("Writing merged raw RNA counts table to: ", output_folder, output_name, "_RNA_raw_counts.txt", sep = ""))
  write.table(rna_counts, file = paste(output_folder, output_name, "_RNA_raw_counts.txt", sep = ""),quote = F, sep = "\t", row.names = T)
  
  #### RNA RPKMs of the fragments 
  # Download gene information using biomaRt 
  ensembl_hg19 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") 
  All_genes_info <- getBM(attributes = c("refseq_mrna", "ensembl_transcript_id"), mart = ensembl_hg19)
  exon_info <- getBM(attributes = c("ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"), mart = ensembl_hg19)
  
  d_transciptID <- merge(d_TruncatedGenes, All_genes_info, by = "refseq_mrna")
  
  all_transcripts <- NULL
  for(transcript in unique(d_transciptID$ensembl_transcript_id)){
    my_transcript <- gtf_trunc[grep(transcript, gtf_trunc$V9),]
    all_transcripts <- rbind(all_transcripts, my_transcript)
  }
  
  ## Get exon ids 
  all_transcripts <- all_transcripts[grep("exon_id", all_transcripts$V9),]
  for(i in 1:nrow(all_transcripts)){
    exon_id <- gsub(".*exon_id ", "", gsub(";", "", all_transcripts$V9[i]))
    ensembl_id <- gsub("gene_id ", "", gsub(";.*", "", all_transcripts$V9[i]))
    all_transcripts$ensembl_exon_id[i] <- exon_id
    all_transcripts$ensembl_fragment_id[i] <- ensembl_id
  }  
  
  ## Get exon sizes
  all_transcripts <- merge(all_transcripts, exon_info, by = "ensembl_exon_id")
  all_transcripts$exon_size <- all_transcripts$exon_chrom_end - all_transcripts$exon_chrom_start
  
  ## Combine exon sizes fragments and complete genes 
  my_agg <- aggregate(all_transcripts$exon_size, by=list(Category=all_transcripts$ensembl_fragment_id), FUN=sum)
  row.names(my_agg) <- my_agg$Category
  my_agg$Category <- NULL
  colnames(my_agg) <- "V2"
  exon_gene_sizes <- read.table(Exon_sizes,sep="\t",header=F,row.names=1)
  my_exon_gene_sizes <- rbind(exon_gene_sizes, my_agg)
  
  raw_read_counts <- rna_counts
  
  ## Calculate RPKMs all samples
  nrsamples <- ncol(raw_read_counts)
  nrrows <- nrow(raw_read_counts)
  
  tab <- matrix(data=NA, nrow=nrrows, ncol=nrsamples)
  
  for (j in 1:nrsamples){
    RPKM = rpkm(raw_read_counts[j], exon_gene_sizes, normalized.lib.sizes=F, log=F)
    for(i in 1:nrrows){
      tab[i,j] = RPKM$V2[i]
    }
  }
  df <- data.frame(tab)
  colnames(df) <- colnames(raw_read_counts)
  rownames(df) <- rownames(raw_read_counts)
  
  write.table(df, file = paste(output_folder, output_name, "_RNA_RPKMs.txt", sep = ""),quote = F, sep = "\t", row.names = T)
  
  print(paste("Merged raw RNA count tables of ", ncol(rna_counts), " samples", sep =""))
  
  ## Create filelist for samples that have truncated genes 
  filelist <- filelist[filelist$Sample %in% colnames(rna_counts),]
}

### DESEQ2
filelist$Subject <- droplevels(as.factor(filelist$Subject))

DE_overview <- data.frame()

print("### Starting differential expression analysis using DESeq2 ###")
print(paste("# Start time =", Sys.time()))

i <- 1

for(subject in levels(filelist$Subject)){
  print(paste("## Calculating differential expression sample ", i, " of ", length(levels(filelist$Subject)), ": ", subject, sep = "" ))
  
  patient_samples <- filelist[filelist$Subject == subject, "Sample"]
  

  colData <- data.frame(condition = ifelse(names(rna_counts) %in% patient_samples, subject, "Control"), row.names = names(rna_counts))
  colData$condition <- factor(colData$condition, levels=c("Control", subject))

  print(colData)
  dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                                colData = colData,
                                design = ~ condition)
  dds <- DESeq(dds, quiet = T)
  res <- results(dds)
  
  write.table(as.data.frame(res), file = paste(DE_output,output_name, "_", subject, "_DE_results.bed", sep = ""), quote = F, sep = "\t", row.names = T)
  
  MA_folder <- paste(DE_output, "MA-plots/", sep = "")
  ifelse(!dir.exists(file.path(MA_folder)), dir.create(file.path(MA_folder), showWarnings = FALSE), FALSE)
  
  pdf(file = paste(MA_folder, output_name, "_", subject,  ".pdf",sep = ""), width = 5, height = 5, pointsize = 10)
  plotMA(res, main=paste("MA-plot_", subject, sep = ""), ylim=c(-4,4),colNonSig = rgb(211,211,211,150, maxColorValue = 255), cex = 0.5)
  dev.off()
  
  # This overview is only used to check the number of DE genes:
  res2 <- res[!is.na(res$padj), ]
  padj <- res2[res2$padj < 0.1,]
  pval <- res2[res2$pvalue < 0.05,]
  log2_neg <- res2[res2$log2FoldChange < -1,]
  log2_pos <- res2[res2$log2FoldChange > 1,]
  
  overview <- data.frame(subject = subject, 
                         replicates = nrow(filelist[filelist$subject == subject,]),
                         total_genes = nrow(res),
                         log2_neg = nrow(log2_neg),
                         log2_pos = nrow(log2_pos),
                         padj_0.1 = nrow(padj),
                         pval_0.05 = nrow(pval))
  DE_overview <- rbind(DE_overview, overview)
  
  i <- i+1
}
print("# Differential expression analysis done")
print(paste("# End time DE =", Sys.time()))
print(paste("# Writing merged Normalized RNA counts table to: ", output_folder, output_name, "_RNA_normalized_counts.txt", sep = ""))
normalized_counts <- counts(dds, normalized = TRUE) 
write.table(normalized_counts, file = paste(output_folder, output_name, "_RNA_normalized_counts.txt", sep = ""),quote = F, sep = "\t", row.names = T)


if(qc == TRUE){
  print("### Starting DESeq2 normalization QC ###")
  
  QC_folder <- paste(DE_output, "QC/", sep = "")
  ifelse(!dir.exists(file.path(QC_folder)), dir.create(file.path(QC_folder), showWarnings = FALSE), FALSE)
  print("Starting QC")
  
  write.table(DE_overview, file = paste(QC_folder, output_name, "_Overview_DE_Genes_sample.txt", sep =""), quote = F, sep = "\t", row.names = F)
  
  # Collect counts only for genes with >1 count in all samples.
  # Many ensembl_gene_ids do not have a count in any of the samples
  rna_counts_subs <- rna_counts[which(rowSums(rna_counts) > 1),]
  normalized_counts_subs <- normalized_counts[which(rowSums(normalized_counts) > 1),]
  
  # Check if normalization worked well (variation between samples should be strongly decreased)
  raw_medians <- apply(rna_counts_subs, 2 , median)
  normalized_medians <- apply(normalized_counts_subs, 2, median)
  
  pdf(file = paste(QC_folder, output_name, "_rna_normalization_QC.pdf", sep = ""), paper = "a4", height = 10)
  
  # Median counts
  par(mar = c(6, 4, 4, 2))
  plot(raw_medians, xlab="sample", ylab="median counts", 
       main = "Median counts before and after normalization (red)", cex = 1.5, col = "blue", axes= F, mgp=c(4.5,1,.5))
  box()
  axis(side = 1, at = 1:length(raw_medians), tck = 1, col.tick = "lightgrey", labels =names(raw_medians), las =2, cex.axis = 0.7)
  axis(side = 2)
  points(normalized_medians, col = "red", cex = 1.5, bg = "red", pch = 21)
  
  # Compare distribution of counts and normalized counts (not essential)
  par(mfrow=c(2,1),oma = c(2,1,1,1))
  boxplot(rna_counts, outline = FALSE, las = 2, 
          ylim = c(0,50),ylab = "Counts per gene",
          main = "Count distribution before normalization", 
          cex.main = 0.8, cex.axis = 0.5, col = "blue")
  boxplot(normalized_counts, outline = FALSE, las = 2, 
          ylim = c(0,50), ylab = "Counts per gene",
          main = "Count distribution after normalization", 
          cex.main = 0.8, cex.axis = 0.5, col = "red")
  par(mfrow=c(1,1))
  
  # Plot the boxplots after removing all genes without counts:
  par(mfrow=c(2,1),oma = c(2,1,1,1))
  boxplot(rna_counts_subs, outline = FALSE, las = 2, 
          ylim = c(0,1000),ylab = "Counts per gene",
          main = "Count distribution before normalization \n (Genes >1 count)", 
          cex.main = 0.8, cex.axis = 0.5, col = "blue")
  boxplot(normalized_counts_subs, outline = FALSE, las = 2, 
          ylim = c(0,1000), ylab = "Counts per gene",
          main = "Count distribution after normalization \n (Genes >1 count)", 
          cex.main = 0.8, cex.axis = 0.5, col = "red")
  par(mfrow=c(1,1))
  
  # Make a barplot of the total counts per sample
  total_counts <- rbind(colSums(rna_counts), colSums(normalized_counts))
  barplot(total_counts/1e6, beside = T, col = c("blue", "red"), ylab = "Total counts (millions)", 
          las = 2, cex.axis = 0.8, cex.main = 0.8, cex.names = 0.5, main = "Total counts per sample \n (Blue = raw, Red = normalized)")
  abline(h = median(total_counts)/1e6)
  
  dev.off()
  
  DE_overview_m <- melt(DE_overview[,c(1,6, 7)])
  pdf(file = paste(QC_folder, output_name, "_DE_Genes_sample.pdf", sep = ""), paper = "a4r", width = 10, height = 8)
  
  DE_results_plot <- ggplot(data = DE_overview_m, aes(x = subject, y = value, fill= variable)) + geom_bar(stat = "identity", position = "dodge") + 
    ggtitle("Number of DE genes per sample") + theme_minimal() +
    theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5), axis.text.x = element_text(angle = 90))
  print(DE_results_plot)
  
  dev.off()
}

#### Clustering
if(plots == TRUE){
  library("pheatmap")
  
  print("### Starting sample clustering ### ")
  
  QC_folder <- paste(DE_output, "QC/", sep = "")
  ifelse(!dir.exists(file.path(QC_folder)), dir.create(file.path(QC_folder), showWarnings = FALSE), FALSE)
  
  
  # Log transformation of dataset, takes a while:
  rld <- rlog(dds, blind=TRUE)
  
  save(rld, file = paste(QC_folder, output_name, "_rlog_counts.RData", sep =""))
  print("## Finished log transformation")
  
  # Heatmap showing sample to sample distance (Euclidean distance)
  
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  df <- data.frame(row.names = row.names(colData(dds)),
                   Patient = gsub(row.names(colData(dds)), pattern = "-RNA-\\d", replacement = ""))
  
  write.table(sampleDistMatrix, file = paste(QC_folder, output_name,"_RNA_sampleDistMatrix.txt", sep =""), sep ="\t")
  
  pdf(paste(QC_folder, output_name, "_RNA_sample_distance_heatmap.pdf", sep =""), pointsize = 12,width = 9, height = 8, onefile=FALSE)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, annotation_col=df, 
           main = "Sample-to-Sample distances")
  dev.off()
  
  # Heatmap gene clustering (Euclidean distance)
  
  variance_per_gene <- rowVars(assay(rld))
  genes_highest_variance <- order(-variance_per_gene)[1:500] 
  my.breaks <- c(seq(-4.4, -1.99, length.out=40),seq(-2, 1.99, length.out=70),seq(2,4.4, length.out=40))
  
  pdf(paste(QC_folder, output_name, "_RNA_sample_clustering_heatmap.pdf", sep =""), onefile=FALSE)
  pheatmap(assay(rld)[genes_highest_variance,], show_rownames=FALSE, annotation_col=df,scale = "row", 
           col = redblue(149), breaks = my.breaks, treeheight_row = 0,
           main = "Sample clustering based on coverage \n 500 genes with highest variance")
  
  
  dev.off()	
  
  # PCA plot
  
  data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))
  
  data$Patient <- gsub(row.names(data), pattern = "-RNA-\\d", replacement = "")
  data$Patient <- factor(data$Patient)
  data$Sample <- gsub(row.names(data), pattern = "\\d\\d\\d\\d\\d-RNA-", replacement = "")
  data$Sample <- factor(data$Sample)
  data$Run <- filelist[,"Run"]
  write.table(data, file = paste(QC_folder, output_name,"_RNA_PCA.txt", sep =""), sep ="\t", quote = F)
  pdf(paste(QC_folder, output_name,"_RNA_PCA.pdf", sep =""), pointsize = 12,width = 9, height = 8)
  
  PCA_plot <- ggplot(data, aes(PC1, PC2, col = Patient, shape = Run)) +
    geom_point(size=3) + geom_text(data = data, aes(x = PC1, y = PC2, label = Patient), nudge_y = 1) +xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    ggtitle("PCA RNA") + 
    theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5))
  
  print(PCA_plot)
  
  dev.off()
  
  print(paste("# End time QC =", Sys.time()))
}

# normalized_counts_m <- melt(normalized_counts)
# graph <- ggplot(data = normalized_counts_m, aes(x = value, fill = Var2)) + geom_density(alpha = 0.3) + coord_cartesian(xlim = c(0, 100)) 
# graph
print("### Differential RNA expression analysis done ###")


