Calculate_TAD_Disruption <- function(Breakpoints = Breakpoints_g,
                                     input_TADs = TAD_file,
                                     cell_types = TADs_Cell_types,
                                     write_output = FALSE,
                                     output_folder = data_folder,
                                     output_name = "Disrupted_TADs_per_patient.txt",
                                     Genes = genes){
  
  print("# Determining which TADs are disrupted")
  TADs <- read.delim(input_TADs, stringsAsFactors = F)
  TADs <- TADs[TADs$TAD_cell %in% cell_types,]
  TADs$TAD_cell <- factor(TADs$TAD_cell)
  
  TADs$TAD_chr <- gsub("chr", "", TADs$TAD_chr)
  TADs$TAD_ID <- paste(TADs$TAD_cell, TADs$TAD_chr, row.names(TADs), sep = "_")
  TADs_G <- GRanges(seqnames = TADs$TAD_chr, IRanges(start = TADs$TAD_start, end = TADs$TAD_end), TAD_cell = TADs$TAD_cell)
  
  disrupted_TADs <- data.frame()
  #number_disrupted_TADs <- data.frame(Cell_type = levels(TADs$TAD_cell))
  
  # First determine which TADs are disrupted in the patient:
  for(Patient in unique(Genes$Patient)){
    #print(Patient)
    
    # Select the breakpoints of the patient
    Breakpoints_Patient <- Breakpoints[which(Breakpoints$Patient == Patient),]
    Breakpoints_Patient_G <- GRanges(seqnames = Breakpoints_Patient$Chr, IRanges(start = Breakpoints_Patient$Breakpoint, end = Breakpoints_Patient$Breakpoint))
    
    # Overlap the breakpoints with all the TADs
    olap_TADs_breakpoints <- findOverlaps(TADs_G,Breakpoints_Patient_G)
    
    # Select the TADs overlapping with breakpoints:
    disrupted_TADs_patient <- TADs[queryHits(olap_TADs_breakpoints),]
    disrupted_TADs_patient$Patient <- Patient
    
    disrupted_TADs_patient$Breakpoint <- start(Breakpoints_Patient_G)[subjectHits(olap_TADs_breakpoints)]
    
    disrupted_TADs_patient$TAD_ID <- factor(disrupted_TADs_patient$TAD_ID)
    disrupted_TADs_patient$TAD_width <- disrupted_TADs_patient$TAD_end - disrupted_TADs_patient$TAD_start
    
    # Next calculate which part of each TAD is disrupted.
    # A TAD that overlaps with a single breakpoint will lead to two "TAD fragments": One before and one after the breakpoint
    # Some TADs overlap with multiple breakpoints, which leads to more than two "TAD Fragments"
    for(TAD in levels(disrupted_TADs_patient$TAD_ID)){
      #print(TAD)
      TAD_Fragments <- disrupted_TADs_patient[disrupted_TADs_patient$TAD_ID == TAD,]
      
      TAD_Fragments$Der_TAD_Start <- TAD_Fragments$TAD_start
      TAD_Fragments$Der_TAD_End <- TAD_Fragments$TAD_end
      
      TAD_Fragments[nrow(TAD_Fragments)+1,] <- TAD_Fragments[nrow(TAD_Fragments),]
      
      for(i in 1:(nrow(TAD_Fragments))){
        #print(i)
        if(i != nrow(TAD_Fragments)){
          # The first TAD fragment starts with the start of the TAD
          # The second TAD fragment will start with the first breakpoint
          TAD_Fragments[i+1, "Der_TAD_Start"] <- TAD_Fragments[i, "Breakpoint"] 
          
          # The first fragment will end with the first breakpoint
          TAD_Fragments[i, "Der_TAD_End"] <- TAD_Fragments[i, "Breakpoint"]
          
          # Each TAD fragment will get a number:
          TAD_Fragments$TAD_Der_ID[i] <- paste(TAD_Fragments$TAD_ID[i], i, sep = "_")
        } else {
          
          # The last fragment will start with the last breakpoint and end with the end of the TAD
          TAD_Fragments[i, "Der_TAD_Start"] <- TAD_Fragments[i, "Breakpoint"]
          TAD_Fragments[i, "Der_TAD_End"] <- TAD_Fragments[i, "TAD_end"]
          TAD_Fragments$TAD_Der_ID[i] <- paste(TAD_Fragments$TAD_ID[i], i, sep = "_")
        }
      }
      # Calculate the width of the TAD fragment:
      TAD_Fragments$Der_Width <- TAD_Fragments$Der_TAD_End - TAD_Fragments$Der_TAD_Start
      
      # TAD lost is the width of the normal/total TAD minus the width of the fragment
      TAD_Fragments$TAD_Lost <- round((TAD_Fragments$TAD_width - TAD_Fragments$Der_Width) / TAD_Fragments$TAD_width * 100,0)
      #print(TAD_Fragments)
      disrupted_TADs <- rbind(disrupted_TADs, TAD_Fragments)
    }
  }
  
  disrupted_TADs_output <- disrupted_TADs[,c(6,5, 1:4, 7:13)]
  
  if(write_output == TRUE){
    write.table(x = disrupted_TADs_output, file = paste(output_folder,output_name, sep = ""), sep = "\t", quote = F, row.names = F)
  }
  return(disrupted_TADs_output)
}

Overlap_Genes_Disrupted_TADs <- function(Disrupted_TADs = TADs_disrupted,
                                         Genes = Genes_SVs,
                                         write_output = FALSE,
                                         output_folder = data_folder,
                                         output_name = "Disrupted_TADs_per_gene.txt"){
  
  # This script determines which genes are located on the previously determined TAD fragments
  
  print("# Determining which genes are located in disrupted TADs")
  Genes_Disrupted_TADs <- data.frame()
  
  for(Patient in unique(Genes$Patient)){
    #print(Patient)
    
    Disrupted_TADs_patient <- Disrupted_TADs[which(Disrupted_TADs$Patient == Patient),]
    
    disrupted_TADs_patient_G <- GRanges(seqnames = Disrupted_TADs_patient$TAD_chr, 
                                        IRanges(start = Disrupted_TADs_patient$Der_TAD_Start, end = Disrupted_TADs_patient$Der_TAD_End))
    Genes_SV_Patient <- Genes[Genes$Patient == Patient,]
    Genes_SV_Patient_G <- GRanges(seqnames = Genes_SV_Patient$chromosome_name, 
                                  IRanges(start = Genes_SV_Patient$transcription_start_site, 
                                          end = Genes_SV_Patient$transcription_start_site+1))
    # Overlap genes with disrupted TADs:
    olap_disrupted_TADs_Genes <- findOverlaps(Genes_SV_Patient_G, disrupted_TADs_patient_G)
    
    Genes_Disrupted_TADs_Patient <- Genes_SV_Patient[queryHits(olap_disrupted_TADs_Genes),]
    Genes_Disrupted_TADs_Patient <- cbind(Genes_Disrupted_TADs_Patient, Disrupted_TADs_patient[subjectHits(olap_disrupted_TADs_Genes),])
    
    # Some genes overlap with multiple (disrupted) TADs in the same cell type.
    Genes_Disrupted_TADs_Patient$Patient_Cell <- paste(Genes_Disrupted_TADs_Patient$ensembl_gene_id, Genes_Disrupted_TADs_Patient$TAD_cell, sep ="_")
    
    # Remove these genes:
    Genes_Disrupted_TADs_cells <- Genes_Disrupted_TADs_Patient[!duplicated(Genes_Disrupted_TADs_Patient$Patient_Cell),]
    
    Genes_Disrupted_TADs_cells$ensembl_gene_id <- as.vector(Genes_Disrupted_TADs_cells$ensembl_gene_id)
    Genes_Disrupted_TADs <- rbind(Genes_Disrupted_TADs, Genes_Disrupted_TADs_cells)
  }
  
  # Only export essential columns:
  Genes_Disrupted_TADs_output <- Genes_Disrupted_TADs[,c("Gene_ID", "Patient", "hgnc_symbol", 
                                                         "chromosome_name","start_position","end_position","strand", "transcription_start_site",
                                                         "TAD_chr","TAD_start","TAD_end","TAD_cell", "TAD_ID","Breakpoint","TAD_width", "Der_TAD_Start", "Der_TAD_End", "TAD_Der_ID",
                                                         "Der_Width", "TAD_Lost", "Patient_Cell")]
  
  #write.table(x = Genes_Disrupted_TADs_cells, file = paste(data_folder, "raw_disrupted_TADs_per_gene.txt", sep = ""), sep = "\t", quote = F)
  if(write_output == TRUE){
    write.table(x = Genes_Disrupted_TADs_output, file = paste(output_folder,output_name, sep = ""), sep = "\t", quote = F, row.names = F)
  }
  return(Genes_Disrupted_TADs)
}

Overview_Disrupted_TADs_per_Gene <- function(Genes_TADs = Genes_Disrupted_TADs,
                                             write_output = FALSE,
                                             output_folder = data_folder,
                                             output_name = "Overview_disrupted_TADs_per_gene.txt"){
  
  # This script generates an overview of the average % of TAD that is lost for each gene (one row per gene)
  
  print("# Making an overview of the disrupted TADs per gene")
  TAD_loss_overview <- data.frame()
  
  for(Patient in unique(Genes_TADs$Patient)){
    #print(Patient)
    # Make an overview for the TAD loss per gene per cell type:
    temp_df <- data.frame(TAD_cell = levels(Genes_TADs$TAD_cell))
    Number_of_Celltypes <-  length(levels(Genes_TADs$TAD_cell))
    Genes_TADs_Patient <- Genes_TADs[Genes_TADs$Patient == Patient,]
    #print(head(Genes_TADs_Patient))
    
    for(gene in levels(factor(Genes_TADs_Patient$hgnc_symbol))){
      #print(gene)
      gene_data <- Genes_TADs_Patient[Genes_TADs_Patient$hgnc_symbol == gene,
                                      names(Genes_TADs_Patient) %in% c("hgnc_symbol", "TAD_cell","TAD_Lost")]
      gene_data <- gene_data[!duplicated(gene_data$TAD_cell),]
      #print(head(gene_data))
      gene_TAD_loss_per_cell <- merge(temp_df, gene_data[,c(2,3)], all.x = T)
      #print(head(gene_TAD_loss_per_cell))
      row.names(gene_TAD_loss_per_cell) <- gene_TAD_loss_per_cell$TAD_cell
      
      # Change all NA's to 0's
      gene_TAD_loss_per_cell$TAD_Lost[is.na(gene_TAD_loss_per_cell$TAD_Lost)] <- 0
      
      # Change the name of the last added column to the gene ID
      names(gene_TAD_loss_per_cell)[names(gene_TAD_loss_per_cell) == "TAD_Lost"] <- paste(Patient, gene, sep ="_")
      
      # Turn the table (genes as rows instead of columns)
      gene_TAD_loss_per_cell_t <- data.frame(t(gene_TAD_loss_per_cell), stringsAsFactors = F)
      gene_TAD_loss_per_cell_t <- gene_TAD_loss_per_cell_t[-1,]
      
      # Calculate in how much % of the cell types the TAD of the gene is disrupted
      gene_TAD_loss_per_cell_t$Cell_Types_TAD_Disrupted <- round(length(which(as.numeric(gene_TAD_loss_per_cell_t) > 0)) / Number_of_Celltypes * 100, 0)
      gene_TAD_loss_per_cell_t$Average_TAD_Loss <- mean(as.numeric(gene_TAD_loss_per_cell_t))
      gene_TAD_loss_per_cell_t$SD_TAD_Loss <- sd(as.numeric(gene_TAD_loss_per_cell_t))
      gene_TAD_loss_per_cell_t$Patient <- Patient
      gene_TAD_loss_per_cell_t$hgnc_symbol <- gene
      gene_TAD_loss_per_cell_t$Gene_ID <- row.names(gene_TAD_loss_per_cell_t)
      
      TAD_loss_overview <- rbind(TAD_loss_overview, gene_TAD_loss_per_cell_t)
      
    }
  }
  
  if(write_output == TRUE){
    write.table(x = TAD_loss_overview, file = paste(output_folder,output_name, sep = ""), sep = "\t", quote = F, row.names = T)
  }
  return(TAD_loss_overview)
}

Calculate_TAD_Loss_Gene <- function(output_folder = data_folder,
                                    output_name = "Overview_disrupted_TADs_per_gene.txt",
                                    overwrite = FALSE,
                                    breakpoints,
                                    genes,
                                    TADs,
                                    cell_types = c("AD", "AO", "BL", "CO", "GM12878", "H1", "HC", "IMR90", "LG", "LI", "LV", "MES","MSC","NPC", "OV", "PA", "PO", "RV", "SB", "SX", "TRO", "CP", "GZ", "DLPFC")){
  
  # breakpoints need to contain columns: Patient - Chr - Breakpoint
  
  output_file = paste(output_folder, output_name, sep = "")
  
  if(file.exists(output_file) == FALSE || overwrite == TRUE ){
    print("Calculating TAD Loss")
    
    TADs_disrupted <- Calculate_TAD_Disruption(Breakpoints = breakpoints,
                                               input_TADs = TADs,
                                               cell_types = cell_types,
                                               write = TRUE,
                                               Genes = genes, 
                                               output_folder = output_folder)
    Genes_Disrupted_TADs <- Overlap_Genes_Disrupted_TADs(Disrupted_TADs = TADs_disrupted,
                                                         Genes = genes,
                                                         write = TRUE, 
                                                         output_folder = output_folder)
    TAD_loss_overview <- Overview_Disrupted_TADs_per_Gene(Genes_TADs = Genes_Disrupted_TADs,
                                                          write = TRUE, output_folder = output_folder)
    
    print("# Finished TAD analysis")
    return(TAD_loss_overview)
    
  } else {
    print("Re-using previously calculated loss per TAD")
    print(paste("Reading: ", output_file))
    
    TAD_loss_overview <- read.delim(output_file)
    return(TAD_loss_overview)
  }
}



### Virtual 4C
Virtual_4C_Disruption <- function(virtual_4c_folder,
                                  output_folder = output_folder,
                                  output_filename = "V4C_Disruption.txt",
                                  Genes = Genes_SVs,
                                  overwrite = TRUE,
                                  re_use = TRUE,
                                  bin_sizes = c(10000,40000)){
  
  # Some HiC matrices have a 10k and some have a 40k resolution (bin size). The V4C profiles (number of columns) for these HiC datasets.
  # This should be the folder structure: [project_folder]/HiC/V4C_2Mb/[resolution]/[chromosome]/[hgnc_symbol].txt
  # The script will loop over each patient, then each gene and then each bin_size
  # Currently only works with V4C files with a certain window around the TSS. These files are much smaller and faster to read.
  output_file <- paste(output_folder, output_filename, sep = "")
  
  if(file.exists(output_file) == FALSE || overwrite == TRUE ){
    
    print(paste("Calculating % disruption of virtual 4C counts for ", nrow(Genes), " genes", sep = ""))
    if(nrow(Genes) > 0){
      V4C_overview <- data.frame()
      
      # counter is just to print progress:
      i <- 1
      warnings <- seq(0.1*nrow(Genes), nrow(Genes), by = nrow(Genes) / 10)
      warning <- 1
      
      for(Patient in unique(Genes$Patient)){
        #print(Patient)
        V4C_disruption_patient <- data.frame()
        for(Gene_ID in Genes$Gene_ID[Genes$Patient == Patient]){
          V4C_loss_gene <- data.frame()
          
          gene_data <- Genes[which(Genes$Gene_ID == Gene_ID),]
          hgnc_symbol <- as.vector(unique(gene_data$hgnc_symbol))
          
          for(resolution in bin_sizes){
            
            virtual_4C_file <- unique(paste(virtual_4c_folder, resolution, "/", gene_data$chromosome_name,"/", hgnc_symbol, ".txt", sep =""))
            
            # Read the pre-computed V4C table for the gene:
            # print(virtual_4C_file)
            if(file.exists(virtual_4C_file) == TRUE){
              virtual_4C_Gene <- read.delim(virtual_4C_file, header = F, stringsAsFactors = F, row.names = 1)
              
              # Usually, when using a V4C file with a specific window around the TSS, the TSS will be in the center column. 
              # However, some genes have less columns if they are close to the start or end of a chromosome.
              bin_gene <- floor(gene_data$transcription_start_site/resolution)
              window <- 2e6
              
              # The first column will be #1 if the gene is close to start of chr
              window_start <- ifelse((bin_gene - window / resolution) > 1, (bin_gene - window / resolution), 1)
              # The final column is the start column plus the total number of columns if the gene is close to the end of the chr
              window_end <- ifelse((bin_gene + window / resolution) < (window_start + ncol(virtual_4C_Gene)), 
                                   (bin_gene + window / resolution), 
                                   (window_start + ncol(virtual_4C_Gene)-1))
              # Column names will contain the bin ID (will be string instead of numeric)
              names(virtual_4C_Gene) <- window_start:window_end
              
              # If the SV starts before the first V4C bin, the fragment start will be set to the number of the first V4C bin              
              bin_fragment_start <- ifelse(floor(gene_data$SV_start / resolution) < window_start, window_start, floor(gene_data$SV_start / resolution))
              # If the SVs ends after the final V4C bin, the fragment end will be set to the number of the final V4C bin
              bin_fragment_end <- ifelse(floor(gene_data$SV_end / resolution) >  window_end, window_end, floor(gene_data$SV_end / resolution))
              
              V4C_counts_gene_window <- rowSums(virtual_4C_Gene[,c(as.character(window_start:window_end))])
              if(bin_fragment_start != bin_fragment_end){
                V4C_counts_fragment <- rowSums(virtual_4C_Gene[,c(as.character(bin_fragment_start:bin_fragment_end))])
                
              } else {
                V4C_counts_fragment <- virtual_4C_Gene[,c(as.character(bin_fragment_start))]
              }
              
              V4C_loss_resolution <- t(data.frame(round((V4C_counts_gene_window - V4C_counts_fragment) / V4C_counts_gene_window * 100, 2)))
              
              row.names(V4C_loss_resolution) <- Gene_ID
              #print(V4C_loss_resolution)
              if(ncol(V4C_loss_gene) > 0){
                V4C_loss_gene <- cbind(V4C_loss_gene, V4C_loss_resolution)
              } else {
                V4C_loss_gene <- V4C_loss_resolution
              }
            }
          }
          
          # some genes may not have V4C data for all cell types and are therefore excluded
          if(ncol(V4C_loss_gene) > 3){
            V4C_disruption_patient <- rbind(V4C_disruption_patient, V4C_loss_gene)
          }
          # counter
          i <- i +1
          if(i > warnings[warning]){
            print(paste(warning*10, "% done", sep = ""))
            warning <- warning+1
          }
        }
        V4C_overview <- rbind(V4C_overview, V4C_disruption_patient)
      }
    }
    dir.create(paste(output_folder, "V4C/", sep = ""), showWarnings = F)
    write.table(file = paste(output_folder,"V4C/", Patient, "_V4C_Disruption.txt", sep = ""), x = V4C_disruption_patient, sep = "\t", quote = F)
  } else {
    print("Re-using previously calculated % of virtual 4C loss")
    print(paste("Reading: ", output_file))
    
    V4C_overview <- read.delim(output_file)
  }
  
  V4C_overview$V4C_Disruption <- rowMeans(V4C_overview)
  V4C_overview$V4C_Disruption_SD <- apply(X = V4C_overview, MARGIN = 1, FUN = "sd")
  return(V4C_overview)
}

  

HiC_Loop_Disruption <- function(loop_file,
                                      genes = genes_SVs,
                                      TSS_window = 5000,
                                      loop_window = 0,
                                      breakpoints){
  
  genes_output <- data.frame(Gene_ID = genes[,"Gene_ID"])
  
  loops <- read.delim(loop_file)
  #head(loops)
  
  breakpoints <- GRanges(seqnames = breakpoints$Chr, IRanges(start = breakpoints$Breakpoint, end = breakpoints$Breakpoint + 1))
  
  TSS <- GRanges(seqnames = genes$chromosome_name, IRanges(start =  genes$transcription_start_site - TSS_window,
                                                           end = genes$transcription_start_site + 1 + TSS_window),
                 ensembl_gene_id = genes$ensembl_gene_id)
  
  
  disrupted_loop_overview <- data.frame()
  all_disrupted_loops <- data.frame()
  
  # First the number of loops are determined for each gene
  
  for(cell_type in unique(loops$cell_type)){
    
    print(cell_type)
    
    # Select only the loops for the specific cell type:
    loops_cell <- loops[which(loops$cell_type == cell_type),]
    #print(head(loops_cell))
    
    fragment_x <- GRanges(seqnames = loops_cell$chr1, IRanges(start = loops_cell$x1-loop_window, end = loops_cell$x2+loop_window))
    fragment_y <- GRanges(seqnames = loops_cell$chr2, IRanges(start = loops_cell$y1-loop_window, end = loops_cell$y2+loop_window))
    
    olap_genes_x <- findOverlaps(fragment_x, TSS)
    olap_genes_y <- findOverlaps(fragment_y, TSS)
    
    genes_x <- cbind(genes[subjectHits(olap_genes_x),], loops_cell[queryHits(olap_genes_x),])
    genes_y <- cbind(genes[subjectHits(olap_genes_y),], loops_cell[queryHits(olap_genes_y),])
    
    genes_with_loops <- rbind(genes_x, genes_y)
    genes_with_loops$hgnc_symbol <- as.vector(genes_with_loops$hgnc_symbol)
    
    genes_with_loops$Gene_Loop <- paste(genes_with_loops$Gene_ID, genes_with_loops$Loop, sep = "_")
    
    # Count the number of loops per gene:
    loops_per_gene <- data.frame(Gene_ID = names(table(genes_with_loops$Gene_ID)), loops = as.numeric(table(genes_with_loops$Gene_ID)))
    names(loops_per_gene)[2] <- paste(cell_type, "_loops", sep ="")
    
    
    
    # Overlap the loops with genes with the breakpoints:
    genes_with_loops_g <- GRanges(seqnames = genes_with_loops$chr1, IRanges(start = genes_with_loops$x2, genes_with_loops$y1))
    
    olap_loops_breakpoints <- findOverlaps(genes_with_loops_g, breakpoints)
    
    disrupted_loops <- cbind(genes_with_loops[queryHits(olap_loops_breakpoints),], as.data.frame(breakpoints[subjectHits(olap_loops_breakpoints),])[-5])
    
    disrupted_loops <- disrupted_loops[!duplicated(disrupted_loops$Gene_Loop),]
    if(nrow(disrupted_loops) > 0 ){
      disrupted_loops_per_gene <- data.frame(Gene_ID = names(table(disrupted_loops$Gene_ID)), Disrupted_loops= as.numeric(table(disrupted_loops$Gene_ID)))
      names(disrupted_loops_per_gene)[2] <- paste(cell_type, "dis_loops", sep = "_")
      
      genes_output <- merge(genes_output, loops_per_gene, by = "Gene_ID", all.x = T)
      genes_output <- merge(genes_output, disrupted_loops_per_gene, by = "Gene_ID", all.x = T)
    }else{
      disrupted_loops_per_gene <- data.frame(Gene_ID = names(table(genes_output$Gene_ID)), Disrupted_loops= NA )
      names(disrupted_loops_per_gene)[2] <- paste(cell_type, "dis_loops", sep = "_")
      
      genes_output <- merge(genes_output, loops_per_gene, by = "Gene_ID", all.x = T)
      genes_output <- merge(genes_output, disrupted_loops_per_gene, by = "Gene_ID", all.x = T)
      
    } 
  }
  
  disrupted_loops <- genes_output[,grep("dis_loops", x = names(genes_output))]
  disrupted_loops <- as.matrix(disrupted_loops)
  disrupted_loops[is.na(disrupted_loops)] <- 0
  genes_output$Total_dis_loops <- rowSums(disrupted_loops)
  
  normal_loops <- genes_output[,-grep("dis_loops", x = names(genes_output))]
  normal_loops <- normal_loops[,grep("loops", x = names(normal_loops))]
  
  normal_loops <- as.matrix(normal_loops)
  normal_loops[is.na(normal_loops)] <- 0
  genes_output$Total_loops <- rowSums(normal_loops)
  
  genes_output$Perc_dis_loops <- genes_output$Total_dis_loops / genes_output$Total_loops * 100
  
  return(genes_output)
  
}


