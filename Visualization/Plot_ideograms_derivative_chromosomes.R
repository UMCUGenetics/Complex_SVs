library(ggbio)
library(GenomicRanges)

## The locations of the centromeres are obtained from the "gaps" table of the UCSC table browser
## http://genome.ucsc.edu/cgi-bin/hgTables > all Tables > gaps

if(exists("gaps") == FALSE){
  gaps <- read.delim(paste(input_folder, "Data/Genome/hg19_gaps.txt", sep = ""), stringsAsFactors = F, check.names = F)
}

centromeres <- gaps[which(gaps$type == "centromere"),]

# This function plots derivative chromosome karyotypes
plot_der_karyotype <- function(Patient,
                               breakpoints = SVs,
                               chr_colors = c(255,0,0,20,
                                              50,205,50,20,
                                              0,0,255,20,
                                              255,140,0,20,
                                              128,0,128,20),
                               xlim = c(0,250e6),
                               chr_height = 0.6,
                               lwd = 0.5,
                               axis = "top"){
  
  print(paste("## Plotting karyotypes derivative chromosomes patient ", Patient, sep = ""))
  
  data(hg19IdeogramCyto, package = "biovizBase")
  
  seqlevels(hg19IdeogramCyto) <- gsub("chr", "", seqlevels(hg19IdeogramCyto))
  hg19 <- hg19IdeogramCyto
  
  breakpoints_patient <- breakpoints[breakpoints$Patient == Patient ,]
  
  # First add colors to each chromosomal fragment (different color for each original chromosome)
  i <- 0
  for(chr in unique(breakpoints_patient$chr)){
    #print(chr)
    #print(i)
    breakpoints_patient$chr_fill[breakpoints_patient$chr == chr] <- rgb(chr_colors[1+i*4], chr_colors[2+i*4], chr_colors[3+i*4], chr_colors[4+i*4], maxColorValue = 255)
    breakpoints_patient$chr_border[breakpoints_patient$chr == chr] <- rgb(chr_colors[1+i*4], chr_colors[2+i*4], chr_colors[3+i*4], 200, maxColorValue = 255)
    i <- i + 1
  }
  
  par(mar=c(0, 4, 0, 0))
  
  plot.new()
  title(Patient, cex.main = 0.8)
  
  deletions <- unique(breakpoints_patient$der_chr)[startsWith("del", x=  unique(breakpoints_patient$der_chr))]
  der_chroms <- unique(breakpoints_patient$der_chr)[order(as.numeric(gsub("der", "", unique(breakpoints_patient$der_chr))), decreasing = T)]
  if(length(deletions) > 0){
    
    der_chroms <- der_chroms[-which(der_chroms %in% deletions)]
  }
  
  ylim <- c(0, ifelse(length(deletions > 0),
                      length(der_chroms)*2+2,
                      length(der_chroms)*2))
  
  ylim[2] <- ifelse(length(der_chroms) < 3, 6, ylim[2])
  
  xlim[2] <- ifelse(max(breakpoints_patient$der_end > 250e6), 300e6,  xlim[2])
  
  plot.window(xlim=c(xlim[1],xlim[2]), ylim= c(ylim[1],ylim[2]))
  
  # Plot the chromosomal axis 
  if(axis == "top"){
    ticks <- seq(xlim[1], xlim[2], by = 50e6)/1e6
    axis(1, at = ticks*1e6, labels = paste(ticks, " Mb", sep = ""), side = 3, cex.axis= 0.8)
  } else if ( axis == "bottom"){
    ticks <- seq(xlim[1], xlim[2], by = 50e6)/1e6
    axis(1, at = ticks*1e6, labels = paste(ticks, " Mb", sep = "") ,side = 1, cex.axis= 0.8)
  }

  y0 <- ifelse(length(deletions) > 0, 2.5, 0.5)
  y0 <- ifelse(length(der_chroms) < 3, 2.5, y0)
  
  for(der_chrom in der_chroms){
    #print(der_chrom)
    
    der_chrom_label <- paste(gsub(pattern = "der", replacement = "der(", x = der_chrom), ")", sep ="")
    mtext(text = der_chrom_label, side = 2, las = 1, at = y0+chr_height/2)
    
    fragments <- breakpoints_patient[which(breakpoints_patient$der_chr == der_chrom),]
    
    number_fragments <- nrow(fragments)
    
    i <- 1
    
    for(der_fragment in fragments$der_fragment[order(fragments$der_fragment)]){
      #print(der_fragment)
      
      fragment_data <- breakpoints_patient[which(breakpoints_patient$der_fragment == der_fragment),]
      names(fragment_data) <- paste("SV", names(fragment_data), sep = "_")
      fragment_g <- GRanges(seqnames = fragment_data$SV_chr, IRanges(start = fragment_data$SV_start, end = fragment_data$SV_end))
      olap_fragment_bands <- findOverlaps(hg19, fragment_g)
      
      bands <- cbind(as.data.frame(hg19[queryHits(olap_fragment_bands)]), as.data.frame(fragment_data[subjectHits(olap_fragment_bands),]))
      
      bands$band_start <- ifelse(bands$SV_type != "Inversion" & bands$SV_type != "Flanking_Inversion",
                                 bands$SV_der_start + (bands$start - bands$SV_start),
                                 bands$SV_der_start + (bands$SV_end - bands$end))
      bands$band_end <- ifelse(bands$SV_type != "Inversion" & bands$SV_type != "Flanking_Inversion",
                               bands$SV_der_start + (bands$end - bands$SV_start),
                               bands$SV_der_start + (bands$SV_end - bands$start))
      bands$band_start <- ifelse(bands$band_start < bands$SV_der_start, bands$SV_der_start, bands$band_start)
      bands$band_end <- ifelse(bands$band_end > bands$SV_der_end, bands$SV_der_end, bands$band_end)
      bands$col <- bands$gieStain
      bands$col  <- factor(bands$col)
      levels(bands$col ) <- list( "grey100" = "gneg",
                                  "brown3" = "stalk", 
                                  "#C1C1C1" = "gpos25",
                                  "#808080" = "gpos50",
                                  "#404040" = "gpos75",
                                  "#000000" = "gpos100",
                                  "brown4" = "acen")
      
      # The centromere bands will be smaller in height:
      bands$y0 <- ifelse(bands$gieStain == "acen", y0+chr_height/5, y0)
      bands$y1 <- ifelse(bands$gieStain == "acen", y0+chr_height-chr_height/5, y0+chr_height)
      
      if(bands$SV_type == "Duplication"){
        Duplicated_bands <- bands
        Duplicated_bands$band_start <- Duplicated_bands$band_start + Duplicated_bands$SV_width
        Duplicated_bands$band_end <- Duplicated_bands$band_end + Duplicated_bands$SV_width
        
        bands <- rbind(bands, Duplicated_bands)
        bands$SV_der_width <- bands$SV_der_width*2
      }
      
      # Plot the bands
      rect(xleft = bands$band_start, xright = bands$band_end, ybottom = bands$y0, ytop = bands$y1, col = as.character(bands$col), border = NA)
      
      outlines <- data.frame(starts = c(min(bands$band_start[bands$band_start < min(bands$band_start[bands$gieStain == "acen" ])]), 
                                        min(bands$band_start[bands$band_start > max(bands$band_start[bands$gieStain == "acen" ])]),
                                        min(bands$band_start[bands$gieStain == "acen" ])),
                             ends = c(max(bands$band_end[bands$band_end <= min(bands$band_start[bands$gieStain == "acen" ])]), 
                                      max(bands$band_end[bands$band_end >= max(bands$band_start[bands$gieStain == "acen" ])]),
                                      max(bands$band_end[bands$gieStain == "acen" ])),
                             ybottom = c(y0,y0, y0+chr_height/5),
                             ytop = c(y0+chr_height, y0+chr_height, y0+chr_height-chr_height/5)
                             )
      
      # Plot the outline and shading of the chromosome
      # rect(xleft = min(bands$band_start), xright = max(bands$band_end), ybottom = y0, ytop = y0+chr_height, 
      #      col = bands$SV_chr_fill, border = bands$SV_chr_border, lwd = 2)
      
      rect(xleft =outlines$starts, xright = outlines$ends, ybottom = outlines$ybottom, ytop = outlines$ytop, 
           col = bands$SV_chr_fill, border = bands$SV_chr_border, lwd = lwd)

      if(i != number_fragments){
        arrows(x0 = max(bands$band_end), y0 = y0-0.4, y1 = y0+chr_height+1,lty = 2, code = 0, col = "darkred", lwd = lwd)
      }
      
      if(max(bands$SV_der_width > 3e6)){
        fragment_start_label <- paste(bands[1,"SV_chr"], as.character(bands[1,"name"]), sep = "")
        
        fragment_end_label <- paste(bands[nrow(bands),"SV_chr"], as.character(bands[nrow(bands),"name"]), sep = "")
        
        text(labels = ifelse(unique(grepl("nver", bands$SV_type)) == TRUE,  fragment_end_label, fragment_start_label), 
             y = y0 + chr_height +0.6, x = min(bands$SV_der_start)+ 2.7e6, col = unique(bands$SV_chr_border), cex = 0.8, srt = 90)
        if(max(bands$SV_der_width > 6e6)){
          text(labels = ifelse(unique(grepl("nver", bands$SV_type)) == TRUE,  fragment_start_label, fragment_end_label), y = y0 + chr_height +0.6, x = max(bands$SV_der_end)- 2.e6, col = unique(bands$SV_chr_border), cex = 0.8, srt = 90)
          
          text(labels = unique(bands$SV_Fragment), 
               y = ifelse(max(bands$SV_der_width > 20e6),  y0 + chr_height + 0.6, y0 - 0.2),
               x = (max(bands$SV_der_end) + min(bands$SV_der_start)) / 2, 
               col = unique(bands$SV_chr_border), cex = 0.8)
        }
        
      }
      i <- i + 1
    }
    y0 <- y0 + 1.8
  }
  
  if(length(deletions) > 0){
    mtext(text = "dels", side = 2, las = 1, at = 0.5 + chr_height/2)
    del <- 0
    for(deletion in deletions){
      
      #print(deletion)
      deletion_data <- breakpoints_patient[which(breakpoints_patient$der_chr == deletion),]
      names(deletion_data) <- paste("SV", names(deletion_data), sep = "_")
      deletion_g <- GRanges(seqnames = deletion_data$SV_chr, IRanges(start = deletion_data$SV_start, end = deletion_data$SV_end))
      olap_deletion_bands <- findOverlaps(hg19, deletion_g)
      
      bands <- cbind(as.data.frame(hg19[queryHits(olap_deletion_bands)]), as.data.frame(deletion_data[subjectHits(olap_deletion_bands),]))
      
      bands$band_start <- ifelse(bands$SV_type != "Inversion" & bands$SV_type != "Flanking_Inversion",
                                 bands$SV_der_start + (bands$start - bands$SV_start),
                                 bands$SV_der_start + (bands$SV_end - bands$end))
      bands$band_end <- ifelse(bands$SV_type != "Inversion" & bands$SV_type != "Flanking_Inversion",
                               bands$SV_der_start + (bands$end - bands$SV_start),
                               bands$SV_der_start + (bands$SV_end - bands$start))
      bands$band_start <- ifelse(bands$band_start < bands$SV_der_start, bands$SV_der_start, bands$band_start)
      bands$band_end <- ifelse(bands$band_end > bands$SV_der_end, bands$SV_der_end, bands$band_end)
      bands$col <- bands$gieStain
      bands$col  <- factor(bands$col)
      levels(bands$col ) <- list( "grey100" = "gneg",
                                  "brown3" = "stalk", 
                                  "#C1C1C1" = "gpos25",
                                  "#808080" = "gpos50",
                                  "#404040" = "gpos75",
                                  "#000000" = "gpos100",
                                  "brown4" = "acen")
      
      
      rect(xleft = bands$band_start + del*25e6, xright = bands$band_end+del*25e6, ybottom = 2.5-2, ytop = 2.5-2+chr_height, col = as.character(bands$col), border = NA)
      
      rect(xleft = min(bands$band_start)+del*25e6, xright = max(bands$band_end)+del*25e6, ybottom = 2.5-2, ytop =2.5-2+chr_height, 
           col = bands$SV_chr_fill, border =bands$SV_chr_border, lwd = lwd)
      
      deletion_label <- ifelse(nrow(bands) > 1,
                               paste(bands$seqnames[1], bands$name[1], "-\n",bands$seqnames[nrow(bands)], bands$name[nrow(bands)] ,sep = ""),
                               paste(bands$seqnames, bands$name,sep = ""))
      
      text(labels =deletion_label,
           x = (min(bands$band_start)+del*25e6 + max(bands$band_end)+del*25e6)/2, y = 2.5-2+chr_height+0.5, cex = 0.8, col = bands$SV_chr_border[1], srt = 90)
      
      del_size <- max(bands$band_end)- min(bands$band_start)
      text(deletion, x = (min(bands$band_start)+del*25e6 + max(bands$band_end)+del*25e6)/2, y = 2.5-2-0.2, cex = 0.8)
      text(labels = paste(as.character(round(del_size/1e6, 3)), " Mb", sep =""), x = (min(bands$band_start)+del*25e6 + max(bands$band_end)+del*25e6)/2, y = 2.5-2-0.4, cex = 0.8)
      
      del <- del+1 
    }
    
  }
  
  
  
}






plot_der_chroms <- function(Patient,
                            breakpoints = SVs,
                            min_size = 1.2,
                            chr_colors = c(255,0,0,
                                           50,205,50,
                                           0,0,255,
                                           255,140,0,
                                           128,0,128),
                            arrow_lwd = 4,
                            arrow_length = 0.15){
  
  print(paste("## Plotting schematic overview derivative chromosomes patient ", Patient, sep = ""))
  
  breakpoints_patient <- breakpoints[breakpoints$Patient == Patient ,]
  breakpoints_patient <- breakpoints_patient[!(breakpoints_patient$type == "Deletion" & breakpoints_patient$width < 10000),] # remove small deletions
  breakpoints_patient_g <- GRanges(seqnames = paste("chr", breakpoints_patient$chr, sep =""), IRanges(start = breakpoints_patient$start, end=breakpoints_patient$end), der_fragment = breakpoints_patient$der_fragment)
  
  # Determine which chromosomal fragments contain a centromere:
  centromeres_g <- GRanges(seqnames = centromeres$chrom, IRanges(start = centromeres$chromStart, end=centromeres$chromEnd))
  breakpoints_patient$centromere <- FALSE
  olap_fragments_centromeres <- findOverlaps(breakpoints_patient_g, centromeres_g)
  breakpoints_patient$centromere[queryHits(olap_fragments_centromeres)] <- TRUE
  
  
  # determine the start and end band name of each fragment
  data(hg19IdeogramCyto, package = "biovizBase")
  hg19 <- hg19IdeogramCyto
  
  olap_chr_bands <- findOverlaps(breakpoints_patient_g, hg19)
  bands_chr <- breakpoints_patient[queryHits(olap_chr_bands),]
  bands_chr$Band <- as.character(hg19$name[subjectHits(olap_chr_bands)])
  
  for(fragment in unique(breakpoints_patient$Fragment)){
    #print(fragment)
    
    bands_fragment <- bands_chr[which(bands_chr$Fragment == fragment),]
    
    breakpoints_patient$band_start[breakpoints_patient$Fragment == fragment] <- bands_fragment[1,"Band"]     
    breakpoints_patient$band_end[breakpoints_patient$Fragment == fragment] <- bands_fragment[nrow(bands_fragment),"Band"]   
  }
  breakpoints_patient$band_start <- paste(breakpoints_patient$chr, breakpoints_patient$band_start, sep ="")
  breakpoints_patient$band_end <- paste(breakpoints_patient$chr, breakpoints_patient$band_end, sep ="")
  
  
  # Each affected chromosome will receive a color listed here:
  i <- 0
  for(chr in unique(breakpoints_patient$chr)){
    #print(chr)
    breakpoints_patient$color[breakpoints_patient$chr == chr] <- rgb(chr_colors[1+i*3], chr_colors[2+i*3], chr_colors[3+i*3], maxColorValue = 255)
    i <- i + 1
  }
  
  # Order the chromosomal fragments based on the derivative chromosome:
  breakpoints_patient <- breakpoints_patient[order(breakpoints_patient$der_fragment),]
  
  # Some patients have many breakpoint junctions on one derivative chromosome which makes it difficult to display.
  # Therefore the der chroms with >6 junctions are split in two "rows"/der chroms in the plot.
  ## Each der chrom can only be displayed in two rows, not more
  breakpoints_patient$split <- 1
  breakpoints_patient2 <- data.frame()
  for(der_chr in unique(breakpoints_patient$der_chr)){
    print(der_chr)
    
    breakpoints_der_chrom <- breakpoints_patient[which(breakpoints_patient$der_chr == der_chr),]
    print(breakpoints_der_chrom)
    if(nrow(breakpoints_der_chrom) > 6){
      split_chr <- data.frame()
      number_fragments <- nrow(breakpoints_der_chrom)
      #print(number_fragments)
      
      split_labels <- c("A", "B", "C", "D", "E", "F")
      
      rows <- ceiling(number_fragments / 6)
      fragment_per_row <- ceiling(number_fragments/rows)
      
      #number_fragments <- number_fragments + rows
      breakpoints_der_chrom$max_split <- rows
      for(i in 1:nrow(breakpoints_der_chrom)){
        print(i)
        
        breakpoints_der_chrom[i,"split"] <- ceiling(i / fragment_per_row)

        if( i %in% seq(from = fragment_per_row, to = number_fragments, by = fragment_per_row) & i != number_fragments){
          split_fragment <- breakpoints_der_chrom[i,]
          split_fragment$split <- split_fragment$split + 1
          split_fragment$der_chr <- paste(split_fragment$der_chr, split_labels[ceiling((i+1) / fragment_per_row)], sep = "_")
        
          breakpoints_der_chrom[i,"der_chr"] <- paste(breakpoints_der_chrom[i,"der_chr"], split_labels[ceiling(i / fragment_per_row)], sep = "_")
          split_chr <- rbind(split_chr, breakpoints_der_chrom[i,])

          split_chr <- rbind(split_chr, split_fragment)
          
        } else {
          breakpoints_der_chrom[i,"der_chr"] <- paste(breakpoints_der_chrom[i,"der_chr"], split_labels[ceiling(i / fragment_per_row)], sep = "_")
          split_chr <- rbind(split_chr, breakpoints_der_chrom[i,])
        }
      }
      breakpoints_patient2 <- rbind(breakpoints_patient2, split_chr)
      
    } else {
      breakpoints_der_chrom$max_split <- 1
      breakpoints_patient2 <- rbind(breakpoints_patient2, breakpoints_der_chrom)
    }
  }
  
  
  breakpoints_patient <- breakpoints_patient2
  breakpoints_patient <- breakpoints_patient[order(breakpoints_patient$der_chr, breakpoints_patient$der_fragment),]
  
  # Each fragment will receive a label based on the length of the fragment.
  # Small fragments will have a size of 1/10 of the der chrom. 
  breakpoints_patient$fragment_length <- ifelse(breakpoints_patient$width < 5e6, "small", "large")
  breakpoints_patient$fragment_length <- ifelse(breakpoints_patient$width > 5e6 & breakpoints_patient$width < 20e6, "medium", breakpoints_patient$fragment_length)
  
  # Deletions will be plotted on a seperate row. If there a deletions the ymax will increase with 1 so the deletions can be plotted on this row
  deletions <- unique(breakpoints_patient$der_chr)[startsWith("del", x=  unique(breakpoints_patient$der_chr))]
  
  plot.new()
  
  ymax <- ifelse(length(deletions) > 0, 
                 (length(unique(breakpoints_patient$der_chr))-length(deletions))*1.5+1,
                 length(unique(breakpoints_patient$der_chr))*1.5)
  ymax <- ifelse(length(unique(breakpoints_patient$der_chr)) < 3, 4.5, ymax)
  
  plot.window(xlim=c(0,10), ylim=c(0,ymax))
  
  title(main = Patient)
  
  chr_number <- 0
  del_number <- 0
  
  # the derivative chromosomes are ordered based on the der chromosome number ("der" is removed and the split chr label (_A _B etc) are removed before ordering).
  # dels automatically become NA and are shown last
  der_chroms <- unique(breakpoints_patient$der_chr)[order(as.numeric(gsub("der|_[^_]+$", "", unique(breakpoints_patient$der_chr))), decreasing = F)]
  
  for(der_chr in der_chroms){
    #print(der_chr)
    
    breakpoints_patient_chr <- breakpoints_patient[which(breakpoints_patient$der_chr == der_chr),]
    #print(breakpoints_patient_chr)
    
    if(der_chr %in% deletions == FALSE){
      # the chroms are plotted from top to bottom:
      y_chr <- ymax - 0.5 - chr_number * 1.5
      #print(y_chr)
      sum_fragments <- 0
      
      breakpoints_patient_chr <- rbind(breakpoints_patient_chr, breakpoints_patient_chr[breakpoints_patient_chr$type == "Duplication",])
      breakpoints_patient_chr <- breakpoints_patient_chr[order(breakpoints_patient_chr$der_fragment),]
      # Minimal length of arrows = min_size (label = "small")
      breakpoints_patient_chr$arrow_length <- ifelse(breakpoints_patient_chr$fragment_length == "small", min_size, 0)
      
      # The remaining length of the der chrom (after subtracting small fragments) is divided over medium and large fragments:
      total_length_long_fragments <- 10 - sum(breakpoints_patient_chr$arrow_length)
      number_large <- length(which(breakpoints_patient_chr$fragment_length == "large"))
      number_medium <- length(which(breakpoints_patient_chr$fragment_length == "medium"))
      
      # Large fragments will be twice the size of medium fragments
      min_size_long_frags <- total_length_long_fragments / (number_large*2 + number_medium)
      breakpoints_patient_chr$arrow_length[which(breakpoints_patient_chr$fragment_length == "large")] <- min_size_long_frags *2
      breakpoints_patient_chr$arrow_length[which(breakpoints_patient_chr$fragment_length == "medium")] <- min_size_long_frags
      
      # Deletions will be plotted with a min_size length
      breakpoints_patient_chr$arrow_length[which(breakpoints_patient_chr$type == "Deletion")] <- min_size
      
      mtext(text = der_chr, side = 2, at = y_chr, las = 1)
      
      #unique(grepl("nversion", bands$SV_type)) == TRUE
      
      for(i in 1:nrow(breakpoints_patient_chr)){
        
        start <- ifelse(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE,
                        sum_fragments,
                        sum_fragments + breakpoints_patient_chr$arrow_length[i])
        
        end <- ifelse(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE,
                      sum_fragments + breakpoints_patient_chr$arrow_length[i],
                      sum_fragments)
        
        #print(i)
        #print(start)
        
        
        # Last fragment 
        if(i == nrow(breakpoints_patient_chr)){
          if(breakpoints_patient_chr$split[i] != max(breakpoints_patient_chr$max_split)){
            if(grepl("nver", breakpoints_patient_chr$type[i]) == TRUE){
              arrows(x0 = start, x1 = end + ((start-end) / 3) * 2,
                     y0 = y_chr, lwd = 2, code = 0,
                     col = as.character(breakpoints_patient_chr$color[i]),
                     lty = 2, length = arrow_length)
              arrows(x0 = end + ((start-end) / 3) * 2, x1 = end,
                     y0 = y_chr, lwd = arrow_lwd, code = 2,
                     col = as.character(breakpoints_patient_chr$color[i]),
                     lty = 1, length = arrow_length)
            } else {
              arrows(x0 = end, x1 = end + ((start-end) / 3) * 2,
                     y0 = y_chr, lwd = 2, code = 0,
                     col = as.character(breakpoints_patient_chr$color[i]),
                     lty = 2, length = arrow_length)
              arrows(x0 = end + ((start-end) / 3) * 2, x1 = start,
                     y0 = y_chr, lwd = arrow_lwd, code = 0,
                     col = as.character(breakpoints_patient_chr$color[i]),
                     lty = 1, length = arrow_length)
            }
          } else {
            arrows(x0 = start, x1 = end,
                   y0 = y_chr, lwd = arrow_lwd, code = 2,
                   col = as.character(breakpoints_patient_chr$color[i]),
                   lty = 1, length = arrow_length)
          }
          
          
          # First fragment
        } else if (i == 1) {
          
          # the der_chrom will start with a sashed line if split == 2
          if(breakpoints_patient_chr$split[i] > 1){
            if(grepl("nver", breakpoints_patient_chr$type[i]) == TRUE){
              
              arrows(x0 = end, x1 = end + ((start-end) / 3) * 2,
                     y0 = y_chr, lwd = 2, code = 0,
                     col = as.character(breakpoints_patient_chr$color[i]),
                     lty = 2, length = arrow_length)
              arrows(x0 = end + ((start-end) / 3) * 2, x1 = start,
                     y0 = y_chr, lwd = arrow_lwd, code = 0,
                     col = as.character(breakpoints_patient_chr$color[i]),
                     lty = 1, length = arrow_length)
            } else {
              arrows(x0 = start, x1 = start + ((end - start) / 3) * 2,
                     y0 = y_chr, lwd = 2, code = 0,
                     col = as.character(breakpoints_patient_chr$color[i]),
                     lty = 2, length =arrow_length)
              arrows(x0 = start + ((end - start) / 3) * 2, x1 = end,
                     y0 = y_chr, lwd = arrow_lwd, code = 2,
                     col = as.character(breakpoints_patient_chr$color[i]),
                     lty = 1, length = arrow_length)
            } 
          }  else {
            arrows(x0 = start, x1 = end,
                   y0 = y_chr, lwd = arrow_lwd, code = 2,
                   col = as.character(breakpoints_patient_chr$color[i]),
                   lty = 1, length = arrow_length)
          }
          
          # else it will start with a normal arrow
        } else {
          arrows(x0 = start, x1 = end,
                 y0 = y_chr, lwd = arrow_lwd, code = 2,
                 col = as.character(breakpoints_patient_chr$color[i]),
                 lty = 1, length = arrow_length)
          
        } 
        
        # BP
        arrows(x0 = ifelse(i > 1, sum_fragments, NA),
               y0 = y_chr-0.5, y1 = y_chr+0.6, lwd = 1.5, code = 0,
               col = "darkred", lty = 2)
        
        # The height of the breakpoint location labels will be the same before and after a junction
        y_label_start <- ifelse(i %% 2 == 0, y_chr-0.35, y_chr-0.2)
        y_label_end <- ifelse(i %% 2 == 0, y_chr-0.2, y_chr-0.35)
        
        if(i != 1){
          if(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE){
            start_label <- format(breakpoints_patient_chr$start[i], big.mark=".", scientific=FALSE)
          } else {
            start_label <- format(breakpoints_patient_chr$end[i], big.mark=".", scientific=FALSE)
          }
        } else if(breakpoints_patient_chr$split > 1){
          start_label <- ""
        } else {
          if(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE){
            start_label <- format(breakpoints_patient_chr$start[i], big.mark=".", scientific=FALSE)
          } else {
            start_label <- format(breakpoints_patient_chr$end[i], big.mark=".", scientific=FALSE)
          }      
        }
        
        
        if(i != nrow(breakpoints_patient_chr)){
          if(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE){
            end_label <- format(breakpoints_patient_chr$end[i], big.mark=".", scientific=FALSE)
          } else {
            end_label <- format(breakpoints_patient_chr$start[i], big.mark=".", scientific=FALSE)
          }
        } else if(breakpoints_patient_chr$split != max(breakpoints_patient_chr$max_split)){
          end_label <- ""
        } else {
          if(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE){
            end_label <- format(breakpoints_patient_chr$end[i], big.mark=".", scientific=FALSE)
          } else {
            end_label <- format(breakpoints_patient_chr$start[i], big.mark=".", scientific=FALSE)
          }      
        }
        
        
        # label start position fragment
        text(x = ifelse(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE, start + 0.5, start - 0.5), 
             labels = ifelse(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE,start_label, end_label),
             y = ifelse(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE, y_label_start,y_label_end), cex = 0.8, font = 3)
        
        #end label
        
        text(x = ifelse(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE, end - 0.5, end + 0.5), 
             labels = ifelse(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE,end_label, start_label),
             y =ifelse(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE, y_label_end,y_label_start), cex = 0.8, font = 3)
        
        ## Band labels
        text(x = ifelse(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE, start + 0.2, start - 0.2), 
             labels = breakpoints_patient_chr$band_start[i],
             y = y_chr+0.35, cex = 0.8, srt = 90, col = breakpoints_patient_chr$color[i])
        
        text(x = ifelse(grepl("nver", breakpoints_patient_chr$type[i]) == FALSE, end - 0.2, end + 0.2), 
             labels = breakpoints_patient_chr$band_end[i],
             y = y_chr+0.35, cex = 0.8, srt = 90, col = breakpoints_patient_chr$color[i])
        
        
        
        # fragment label
        text(x = (start + end) / 2, y = y_chr + 0.5, labels = breakpoints_patient_chr$Fragment[i], font = 2)
        
        # SV type label
        text(x = (start + end) / 2, y = y_chr - 0.5, labels = breakpoints_patient_chr$type[i], font = 3, cex = 0.8)
        
        # Telomere start
        if(i == 1){
          if(breakpoints_patient_chr$type[i] != "Deletion"){
            if(breakpoints_patient_chr$split[i] == 1){
              rect(xleft = min(c(start,end)) - 0.15, xright = min(c(start,end)), ybottom = y_chr-0.3, ytop = y_chr+0.3, col =  as.character(breakpoints_patient_chr$color[i]), border = NA)
            } 
          }
        }
        
        # Telomere end
        if(i == nrow(breakpoints_patient_chr)){
          if(breakpoints_patient_chr$type[i] != "Deletion"){
            if(breakpoints_patient_chr$split[i] == max(breakpoints_patient_chr$max_split)){
              rect(xleft = max(c(start,end)), xright = max(c(start,end))+0.15, ybottom =y_chr-0.3, ytop = y_chr+0.3, col = as.character(breakpoints_patient_chr$color[i]), border = NA)
            }
          }
        }
        
        # centromere
        if(breakpoints_patient_chr$centromere[i] == TRUE){
          points(x = (start+end)/2, y = y_chr, cex = 4, col = as.character(breakpoints_patient_chr$color[i]), pch = 19)
        }
        sum_fragments <- sum_fragments + breakpoints_patient_chr$arrow_length[i]
      }
      
      # Increase the chr_number so the next chromosome will be plotted below this chromosome: 
      chr_number <- chr_number+1
      
    } else {
      
      y_del <- ymax - 0.5 - chr_number * 1.5
      
      breakpoints_patient_chr$arrow_length[which(breakpoints_patient_chr$type == "Deletion")] <- min_size
      
      
      arrows(x0 = del_number*breakpoints_patient_chr$arrow_length, 
             x1 = del_number*breakpoints_patient_chr$arrow_length+breakpoints_patient_chr$arrow_length,
             y0 = y_del, lwd = 2, code = 2,
             col = as.character(breakpoints_patient_chr$color),
             lty = 6, length = arrow_length)
      
      arrows(x0 = del_number*breakpoints_patient_chr$arrow_length,
             y0 = y_del-0.1, y1 = y_del+0.1, lwd = 2, code = 0,
             col = as.character(breakpoints_patient_chr$color),
             lty = 1, length = arrow_length)
      
      
      
      del_length <- breakpoints_patient_chr$width
      
      # name deletion
      text(x = (del_number*breakpoints_patient_chr$arrow_length + del_number*breakpoints_patient_chr$arrow_length+breakpoints_patient_chr$arrow_length) / 2, 
           y = y_del + 0.3, labels = breakpoints_patient_chr$der_chr, font = 2, cex =0.9)
      
      text(x = (del_number*breakpoints_patient_chr$arrow_length + del_number*breakpoints_patient_chr$arrow_length+breakpoints_patient_chr$arrow_length) / 2, 
           y = y_del + 0.1, labels = breakpoints_patient_chr$band_start, font = 3, cex = 0.8)
      
      
      # size deletion
      text(x = (del_number*breakpoints_patient_chr$arrow_length + del_number*breakpoints_patient_chr$arrow_length+breakpoints_patient_chr$arrow_length) / 2, 
           y = y_del - 0.5, labels =  paste("(",as.character(round(breakpoints_patient_chr$width / 1e6, 3)), "Mb)", sep = ""), font = 1, cex = 0.8)
      
      # start-end deletion
      text(x = del_number*breakpoints_patient_chr$arrow_length+0.35, 
           y = y_del - 0.2, labels = format(breakpoints_patient_chr$start, big.mark=".", scientific=FALSE), font = 3, cex = 0.8)
      
      # start-end deletion
      text(x = del_number*breakpoints_patient_chr$arrow_length+breakpoints_patient_chr$arrow_length-0.35, 
           y = y_del - 0.35, labels = format(breakpoints_patient_chr$end, big.mark=".", scientific=FALSE), font = 3, cex = 0.8)
      
      
      del_number <- del_number + 2 
      
      
    }
  }
}

