# Rscript filter_freec.R [FILELIST] [OUTPUT_FOLDER]

# this script uses a filelist to filter multiple CNV files against each other. 
# The first two columns in the filelist should be: sample	freec_file
# the other columns are used for filtering. Each sample will be filtered against the files in the columns on the same row (column names in filelist are used in the output file).
# the sample column should contain unique sample names (sample name will be used for output file)

# example filelist:

# sample	freec_file	father	mother	population
# MP01-3554-DNA-1	[DIR]/freec/MP01-3554-DNA-1_dedup.realigned.bam_CNVs	[DIR]/freec/MP01-1961-DNA-1_dedup.realigned.bam_CNVs	[DIR]/freec/MP01-2428-DNA-1_dedup.realigned.bam_CNVs

# the script will determine overlap between CNV files based on a minimum % of overlap 
# it filters gains and losses separately (a gain overlapping with a loss is not filtered for example) and merges them at the end.
# (!) therefore input files need to contain a column (named "label", "type" or column #5) with "gain" or "loss"

library(GenomicRanges)

args <- commandArgs(trailingOnly=TRUE)
freec_filelist <- args[1]
output_folder <- args[2]

ifelse(!dir.exists(file.path(output_folder)), dir.create(file.path(output_folder), showWarnings = TRUE), FALSE)

findOverlapping <- function(x, y, min_overlap = 0.3, max_overlap = 0.9, max_size = 1e6) {
  # x and y should be GRanges objects
  # will first find overlap between x and y and subsequently determine the number of overlapping bps per hit.
  # number of overlapping bp will be divided by the length of one of the overlapping fragments (longest)
  
  hits <- findOverlaps(x, y)
  xhits <- x[queryHits(hits)]
  yhits <- y[subjectHits(hits)]
  
  # determine how much overlap there is between the two ranges:
  frac <- width(pintersect(xhits, yhits)) / pmax(width(xhits), width(yhits))
  
  # determine the minimal overlap that is necessary to merge the two ranges.
  # depends on the size of the largest fragment (between min_overlap and max_overlap)
  minfrac <- min_overlap + floor(pmax(width(xhits), width(yhits)) / max_size)
  minfrac[which(minfrac > max_overlap)] <- max_overlap
  
  merge <- frac >= minfrac
  return(xhits[merge])
}

flag_FREEC <- function(x, y, name, chr = 1, start = 2,end = 3){
  # x and y need to be dataframes (such as FREEC output)
  # function is used to run findOverlapping function. 
  # add a column named [name] with 1's if overlapping and 0's if there is no overlap between x and y
  
  x$ID <- paste(x[,chr], x[,start], x[,end], sep = "_")
  x_g <- GRanges(seqnames = x[,chr], IRanges(start = x[,start], end = x[,end]), ID = x$ID)
  y_g <- GRanges(seqnames = y[,chr], IRanges(start = y[,start], end = y[,end]))
  
  olap <- findOverlapping(x_g, y_g)
  
  #print(olap)
  x[,name] <- 0
  x[,name] <- ifelse(x$ID %in% olap$ID, 1, 0)
  
  return(x)
}


# read the filelist necessary for filtering:
freec_files <- read.delim(freec_filelist, header = T, stringsAsFactors = F)
print(head(freec_files))

for( i in 1:nrow(freec_files)){
  # loop over each sample, load sample data and subsequently loop over the columns in the filter filelist
  sample <-  freec_files[i,"sample"]
  print(sample)
  freec_file <- freec_files[i,"freec_file"]
  
  input_data <- read.delim(freec_file, header = F, stringsAsFactors = F)
  
  input_gains <- input_data[input_data[,5] == "gain",]
  input_losses <- input_data[input_data[,5] == "loss",]
  
  for(filter_name in names(freec_files)[3:ncol(freec_files)]){
    print(filter_name)
    
    filter_data <- read.delim(freec_files[i,filter_name], header = F)
    if(is.na(as.numeric(gsub("chr", replacement = "", x = "chrX"))) == TRUE){
      filter_data <- read.delim(freec_files[i,filter_name], header = T)
    }
    
    # select the column contain the gain or loss values. This column can be named "type" or "label" or it can be column number 5.
    if(length(grep("type", x = names(filter_data))) == 1){
      gain_loss_column <- "type"
    } else if(length(grep("label", x = names(filter_data))) == 1) {
      gain_loss_column <- "label"
    } else {
      gain_loss_column <- 5
    }
    
    gains_filter <- filter_data[which(filter_data[,gain_loss_column] == "gain"),]
    losses_filter <- filter_data[which(filter_data[,gain_loss_column] == "loss"),]
    
    input_gains <- flag_FREEC(x = input_gains, y = gains_filter,name = filter_name)
    input_losses <- flag_FREEC(x = input_losses, y = losses_filter,name = filter_name)
  }
  
  # merge the flagged gains and losses
  output <- rbind(input_gains, input_losses)
  
  # reorder the output file based on chr - start - end 
  output$V1 <- factor(output$V1, levels = c(1:22, "X", "Y"))
  output <- output[order(output$V1, output$V2, output$V3),]
  
  # remove the rows that have a hit in at least one of the filter columns
  # first 6 columns are the standard freec columns + added ID column
  if(ncol(output) > 7){
    # if the input file is only filtered against multiple files:
    output_filtered <- output[which(rowSums(output[,7:ncol(output)]) < 1),]
  } else {
    # if the input file is only filtered against 1 file:
    output_filtered <- output[which(output[,7] < 1),]
  }
  
  # export the raw output and the filtered output
  write.table(file = paste(output_folder, sample, "_freec_flagged.txt", sep = ""), x = output, sep = "\t", quote = F, row.names = F)
  write.table(file = paste(output_folder, sample, "_freec_filtered.txt", sep = ""), x = output_filtered, sep = "\t", quote = F, row.names = F)
}

