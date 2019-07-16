# This function determines which PCHiC loops overlap with breakpoints
PCHiC_disruption <- function(Input_genes = Genes,
                             Genes_near_SVs = Genes_SVs,
                             PCHiC_data_folder = PCHiC_Folder,
                             cell_types = PCHiC_Celltypes,
                             Input_Breakpoints = Breakpoints,
                             PCHiC_output_folder = paste(output_folder,"PCHiC/", sep = ""),
                             PCHiC_output_filename = "PCHiC_Disruption.txt",
                             overwrite = TRUE){
  
  print("Calculating the % of PCHiC links that are disrupted by SVs")
  dir.create(PCHiC_output_folder, showWarnings = F)
  
  output_file <- paste(PCHiC_output_folder, PCHiC_output_filename, sep = "")
  
  PCHiC_overview <- data.frame()
  
  if(file.exists(output_file) == FALSE || overwrite == TRUE ){
    
    Input_genes <- Input_genes[,c("hgnc_symbol","ensembl_gene_id", "chromosome_name","start_position","end_position","transcription_start_site")]
    
    output <- data.frame()
    
    # First determine which gene is located within the PCHiC bait by overlapping the baits with the transcription start sites:
    Genes_G <- GRanges(seqnames = Input_genes$chromosome_name,
                       IRanges(start = Input_genes$transcription_start_site,
                               end = Input_genes$transcription_start_site + 1),
                       ensembl_gene_id = Input_genes$ensembl_gene_id)
    
    # The script loops over each cell type and reads the PCHiC Data for a cell type. 
    # Then it calculates how many interactions each gene has in that cell type (overlap baits - genes). 
    # Subsequently it selects the genes near the SVs and determines how many interactions are disrupted by the breakpoints (overlap loops - breakpoints)
    for(cell_type in cell_types){
      
      #print(cell_type)
      
      PCHiC_file <- paste(PCHiC_data_folder, "PCHiC_", cell_type, ".txt", sep = "")
      
      if(file.exists(PCHiC_file)){
        print(paste("# Reading ", PCHiC_file, sep = ""))
        PCHiC_data <- read.delim(PCHiC_file, stringsAsFactors = F)
        PCHiC_data$bait_chr <- gsub(pattern = "chr", replacement = "", x = PCHiC_data$bait_chr)
        PCHiC_data$PIR_chr <- gsub(pattern = "chr", replacement = "", x = PCHiC_data$PIR_chr)
        
        PCHIC_Baits_G <- GRanges(seqnames = gsub(pattern = "chr", replacement = "", x = PCHiC_data$bait_chr),
                                 IRanges(start = PCHiC_data$bait_start, 
                                         end = PCHiC_data$bait_end), 
                                 bait_ID =  PCHiC_data$bait_ID, PIR_ID = PCHiC_data$PIR_ID)
        
        overlap_PCHIC_Genes <- findOverlaps(PCHIC_Baits_G, Genes_G)
        
        PCHiC_Genes <- PCHiC_data[queryHits(overlap_PCHIC_Genes),]
        PCHiC_Genes <- cbind(PCHiC_Genes, Input_genes[subjectHits(overlap_PCHIC_Genes),])
        
        # only select the cis-interactions (bait and capture on same chromosome)
        PCHiC_Genes <- PCHiC_Genes[which(PCHiC_Genes$bait_chr == PCHiC_Genes$PIR_chr),]
       
        # Count the number of filtered PCHiC interactions for each gene:
        number_interactions_gene <- data.frame(ensembl_gene_id = names(table(PCHiC_Genes$ensembl_gene_id)),
                                               PCHiC_Interactions = as.numeric(table(PCHiC_Genes$ensembl_gene_id)))
        
        Genes_PCHiC <- merge(Input_genes, number_interactions_gene, by = "ensembl_gene_id", all.x = T)
        #print(head(Genes_PCHiC))
        
        # Only select the PCHiC data for genes near the SVs:
        Genes_subset <-  Genes_near_SVs[,c("ensembl_gene_id", "Gene_ID")]
        Genes_PCHiC_Patient <- merge(Genes_subset, number_interactions_gene, by = "ensembl_gene_id", all.x = T)
        PCHiC_data <- PCHiC_Genes[which(PCHiC_Genes$ensembl_gene_id %in% as.vector(Genes_PCHiC_Patient$ensembl_gene_id)),]
        
        # If the PIR is upstream of the gene/bait: start of the loop = end of PIR & end of the loop = start of the bait
        # If the PIR is downstream of the gene/bait: start of the loop = end of the bait & end of the loop = start of the PIR
        Loops_G <- GRanges(seqnames = PCHiC_data$bait_chr,
                           IRanges(start = ifelse(PCHiC_data$bait_start < PCHiC_data$PIR_start, PCHiC_data$bait_end, PCHiC_data$PIR_end), 
                                   end = ifelse(PCHiC_data$bait_start < PCHiC_data$PIR_start, PCHiC_data$PIR_start, PCHiC_data$bait_start)),
                           bait_ID =  PCHiC_data$bait_ID, PIR_ID = PCHiC_data$PIR_ID)
        
        Breakpoints_g <- GRanges(seqnames = Input_Breakpoints$Chr, IRanges(start=Input_Breakpoints$Breakpoint, 
                                                                     end = Input_Breakpoints$Breakpoint+1),
                                 Patient = Input_Breakpoints$Patient)
        
        # Overlap the PCHiC loops with the breakpoints:
        olap_loops_breakpoints <- findOverlaps(Loops_G, Breakpoints_g)
        disrupted_PCHiC_loops <- PCHiC_data[queryHits(olap_loops_breakpoints),]
        disrupted_PCHiC_loops <- cbind(disrupted_PCHiC_loops, Input_Breakpoints[subjectHits(olap_loops_breakpoints),])
        
        disrupted_PCHiC_loops$Gene_ID <- paste(disrupted_PCHiC_loops$Patient, disrupted_PCHiC_loops$hgnc_symbol, sep = "_")
        
        # Some interactions overlap with multiple breakpoints and therefore appear more than once in list. Remove duplicates:
        disrupted_PCHiC_loops$Loop_ID <- paste(disrupted_PCHiC_loops$Gene_ID, disrupted_PCHiC_loops$ID, sep = "_")
        disrupted_PCHiC_loops <- disrupted_PCHiC_loops[!duplicated(disrupted_PCHiC_loops$Loop_ID),]
        disrupted_PCHiC_loops$ensembl_gene_id <- as.vector(disrupted_PCHiC_loops$ensembl_gene_id)
        
        disrupted_PCHIC_Gene <- data.frame(Gene_ID = names(table(disrupted_PCHiC_loops$Gene_ID)),
                                           Disrupted_PCHiC_interactions = as.numeric(table(disrupted_PCHiC_loops$Gene_ID)))
        
        # Add the disrupted PCHiC data for this cell type to the genelist:
        disrupted_PCHiC_cell <- merge(Genes_PCHiC_Patient, disrupted_PCHIC_Gene, by = "Gene_ID", all.x =T)
        disrupted_PCHiC_cell <- disrupted_PCHiC_cell[,c("Gene_ID", "PCHiC_Interactions", "Disrupted_PCHiC_interactions")]
        
        # Change the names of the columns to the specific cell type:
        names(disrupted_PCHiC_cell) <- c("Gene_ID", paste(cell_type, "_Interactions", sep = ""),
                                         paste(cell_type, "_Disrupted", sep = ""))
        
        disrupted_PCHiC_cell <- disrupted_PCHiC_cell[!is.na(disrupted_PCHiC_cell[,2]),] 
        
        write.table(x = disrupted_PCHiC_cell, file = paste(PCHiC_output_folder, cell_type, "_PCHiC_Disruption.txt"), sep = "\t", quote = F, row.names = F)
        
        
        
        # # Add the cell type PCHiC data to the output overview table:
        # if(nrow(PCHiC_overview) < 1){
        #   PCHiC_overview <- disrupted_PCHiC_cell
        # } else {
        #   PCHiC_overview <- merge(PCHiC_overview, disrupted_PCHiC_cell, by = "Gene_ID", all.x = T)
        # }
      }
      #print(dim(PCHiC_overview))
      
    } 
    
    for(cell_type in cell_types){
      PCHiC_cell <- read.delim(paste(PCHiC_output_folder, cell_type, "_PCHiC_Disruption.txt"), stringsAsFactors = F)

      if(nrow(PCHiC_overview) < 1){
        PCHiC_overview <- PCHiC_cell
      } else {
        PCHiC_overview <- merge(PCHiC_overview, PCHiC_cell, by = "Gene_ID", all = T)
      }
      #print(head(PCHiC_overview))
      
    }
    PCHiC_overview[is.na(PCHiC_overview)] <- 0
    PCHiC_overview$PCHiC_Interactions <- rowSums(PCHiC_overview[,grep("_Interactions", names(PCHiC_overview))])
    PCHiC_overview$PCHiC_Disruptions <- rowSums(PCHiC_overview[,grep("_Disrupted", names(PCHiC_overview))])
    write.table(x = PCHiC_overview, file = output_file, sep = "\t", quote = F, row.names = F)
    
    # # Change all NA's to 0's
    # PCHiC_overview[is.na(PCHiC_overview)] <- 0
    # 
    # # Count the total number of PCHiC interactions and the total number of disrupted PCHiC interactions per gene:
    # PCHiC_overview$PCHiC_Interactions <- rowSums(PCHiC_overview[,grep("_Interactions", names(PCHiC_overview))])
    # PCHiC_overview$PCHiC_Disruptions <- rowSums(PCHiC_overview[,grep("_Disrupted", names(PCHiC_overview))])
    # 
    # write.table(x = PCHiC_overview, file = output_file, sep = "\t", quote = F, row.names = F)
    
  } else {
    print("Re-using previously calculated % of PCHiC links disrupted by SVs")
    print(paste("Reading: ", output_file))
    
    PCHiC_overview <- read.delim(output_file)
  }
  
  return(PCHiC_overview)
}








  
