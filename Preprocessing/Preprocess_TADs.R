# This script is used to process raw files containing TAD boundaries or TAD coordinates to one file containing the TAD boundaries of all cell types
library(GenomicRanges)
library(ggplot2)
library(ggthemes)

options(scipen=999)

## ADD PATH TO OUTPUT AND PROJECT FOLDER
output_folder <- ""
dir.create(output_folder, showWarnings = F)
project_folder <- ""

# txt file containing the lenghts of the chromosomes (hg19). Column 1 = chr ID (chr[.]) Column 2 = length (bp). Can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes 
chromosome_sizes <- read.delim(paste(project_folder, "Common_data/Genes/hg19_chr_sizes.txt", sep = ""), header = F)

# The processed TADs will be stored in Output_TADs
Output_TADs <- data.frame()


## Schmitt et al provided a Excel file containing the TAD boundaries of 21 cell types. Can be downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87112.
# The boundaries of each cell type are stored in separate sheets. 
# Each sheet should be copied to a tab separated txt file with prefix "schmitt2016_TADS_"[celltype].txt . These txt files should be stored in the input_folder_schmitt
# Should contain the TAD boundaries in txt format (chr - start- end)
input_folder_schmitt <- paste(project_folder, "Common_data/TADs/Raw/Schmitt_2016/", sep = "")

# Rewrite the TAD boundaries from Schmitt et al to a txt file containing TADs (Chr - start - end)
TADs <- data.frame()
TAD_Boundaries <- data.frame()

TAD_files <- list.files(path = input_folder_schmitt, pattern = "schmitt2016_TADS_",full.names = T)

for(file in TAD_files){
  
  TADs <- data.frame()
  
  print(file)
  
  # Retrieve the cell type from the file name:
  cell_type <- tail(unlist(strsplit(file, split = "_")),1)
  cell_type <- gsub(pattern = ".txt", replacement = "", x = cell_type)
  
  boundaries <- read.delim(file, header = F, stringsAsFactors = F)

  boundaries$V1 <- gsub(pattern = "chr", x = boundaries$V1, replacement = "")
  names(boundaries) <- c("chr", "start", "end")
  
  TAD_Boundaries_cell <- boundaries
  TAD_Boundaries_cell$Cell_type <- cell_type
  TAD_Boundaries <- rbind(TAD_Boundaries, TAD_Boundaries_cell)
  
  # The boundaries have a width of 40kb. Select the center of each boundary:
  # boundaries$boundary <- (boundaries$start + boundaries$end) / 2
  # boundaries$chr <- factor(boundaries$chr)
  
  for(chr in unique(boundaries$chr)){
    # Select the size of the chromosome (for the last TAD of the chromosome)
    chr_size <- chromosome_sizes[chromosome_sizes$V1 == paste("chr",chr, sep = ""), "V2"]
    
    boundaries_chr <- boundaries[boundaries$chr == chr,]
    
    TAD_Starts <- c(1, boundaries_chr$end)
    TAD_Ends <- c(boundaries_chr$start, chr_size)
    
    TADs_chr <- data.frame(TAD_chr = chr, TAD_start = TAD_Starts, TAD_end = TAD_Ends, TAD_cell = cell_type)
    
    TADs <- rbind(TADs, TADs_chr)
  }
  Output_TADs <- rbind(Output_TADs, TADs)
  
  write.table(paste(output_folder, "TADs_",cell_type,".txt", sep = ""), x = TADs, sep = "\t", quote = F, row.names = F, col.names = F)
}


## Won et al used HiC to determine TADs in two brain cell types. 
# The TAD files (GSE77565_GZ_TAD.bed.gz	and GSE77565_CP_TAD.bed.gz) can be downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77565
# Headers need to be adjusted before merging to Schmitt TADs
Won_2016_CP <- read.delim(paste(project_folder, "Common_data/TADs/Raw/Won_2016/GSE77565_CP_TAD.bed",sep = ""), header = F)
Won_2016_CP$TAD_cell <- "CP"
names(Won_2016_CP) <- names(Output_TADs)
Won_2016_GZ <- read.delim(paste(project_folder, "Common_data/TADs/Raw/Won_2016/GSE77565_GZ_TAD.bed", sep = ""), header = F)
Won_2016_GZ$TAD_cell <- "GZ"
names(Won_2016_GZ) <- names(Output_TADs)

Won_TADs <- rbind(Won_2016_CP, Won_2016_GZ)
Won_TADs$TAD_chr <- gsub(pattern = "chr", replacement = "", x = Won_TADs$TAD_chr)
# This will merge the TAD boundaries (not the TADs themselves)
TAD_Boundaries_all <- TAD_Boundaries
for(chr in levels(Won_TADs$TAD_chr)){
  print(chr)
  data_chr <- Won_TADs[Won_TADs$TAD_chr == chr,]
  for(cell_type in levels(factor(data_chr$TAD_cell))){
    data <- data_chr[data_chr$TAD_cell == cell_type,]
    boundaries <- data.frame(chr = gsub(pattern = "chr", replacement = "", x = chr), 
                             start = c(data$TAD_start[1], data$TAD_end[-length(data$TAD_end)]), 
                             end = c(data$TAD_start[1]+1,data$TAD_start[-1]), Cell_type = cell_type)
    
    TAD_Boundaries_all <- rbind(TAD_Boundaries_all, boundaries)
  }
}
# Add the TADs from Won et al to the other TADs
Output_TADs <- rbind(Output_TADs,Won_TADs)


## Javierre et al determined the TADs of 17 blood cell types. 
# Can be downloaded from: https://osf.io/pg3na/ (extract file TAD_definitions.tar.gz to the folder specified below)
TAD_files_Javierre_2016 <- list.files(path = paste(project_folder, "Common_data/TADs/Raw/Javierre_2016/", sep = ""), pattern = "mean_merged.bed",full.names = T)
for(file in TAD_files_Javierre_2016){
  print(file)
  
  # Retrieve the cell type from the file name:
  filename <- tail(unlist(strsplit(file, split = "/")),1)
  cell_type <- unlist(strsplit(filename, split = "_"))[2]

  TADs_cell <- read.delim(file, header = T)
  
  # Remove the last column
  TADs_cell <- TADs_cell[, -4]
  TADs_cell$TAD_cell <- cell_type
  names(TADs_cell) <- names(Output_TADs)
  Output_TADs <- rbind(Output_TADs, TADs_cell)
}


## Wang et al determined the TADs in adult dorsolateral prefrontal cortex (DLPFC)
# TADs ("DER-18_TAD_adultbrain.bed)" can be downloaded from: http://resource.psychencode.org/

TADs_Wang_DLPFC <- read.delim(paste(project_folder, "Common_data/TADs/Raw/Wang_2018/DER-18_TAD_adultbrain.bed",sep = ""), header = F)
TADs_Wang_DLPFC$TAD_cell <- "DLPFC"
names(TADs_Wang_DLPFC) <- names(TADs)
TADs_Wang_DLPFC$TAD_chr <- gsub(x = TADs_Wang_DLPFC$TAD_chr, pattern = "chr", replacement = "")

# Add the TADs from Wang et al to the other TADs
Output_TADs <- rbind(Output_TADs,TADs_Wang_DLPFC)

summary(factor(Output_TADs$TAD_cell))

write.table(paste(output_folder, "TADs.txt", sep = ""), x = Output_TADs, sep = "\t", quote = F)

#write.table(paste(output_folder, "schmitt2016_TADs.txt", sep = ""), x = TADs, sep = "\t", quote = F)
#write.table(paste(output_folder, "schmitt2016_Boundaries.txt", sep = ""), x = TAD_Boundaries, sep = "\t", quote = F)
#write.table(paste(output_folder, "Boundaries.txt", sep = ""), x = TAD_Boundaries_all, sep = "\t", quote = F)

pdf(paste(project_folder,"Common_data/TADs/TADs_per_Cell.pdf", sep = ""), width = 8, height = 6)
ggplot(data = Output_TADs, aes(x = TAD_cell)) + geom_histogram(stat="count", col = "black", fill = "gray", width = 0.8) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Cell Type", y = "Number of TADs") + ggtitle("Number of TADs per cell type")
dev.off()
