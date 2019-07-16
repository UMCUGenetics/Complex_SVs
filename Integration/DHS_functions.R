DHS_disruption <- function(Genes_near_SVs = Genes_SVs,
                           Input_Breakpoints = Breakpoints,
                           overwrite = TRUE,
                           output_folder = "",
                           DHS_dataset){
  
  if(output_folder != ""){
    # The disrupted connections will be written (not the number of disrupted interactions per gene)
    output_filename <- "Disrupted_DHS_connections"
    output_file_connections <- paste(output_folder, output_filename, ".txt", sep = "")
    output_file_genes <- paste(output_folder, "Disrupted_DHS_connection_counts.txt", sep = "")
  }
  
  print("## Calculating disruptions of DHS-gene connections")
  Genes_near_SVs <- Genes_near_SVs[,c("Gene_ID","hgnc_symbol","ensembl_gene_id","chromosome_name","start_position","end_position","transcription_start_site")]
  
  if(file.exists(output_file) == FALSE || overwrite == TRUE ){
  
    print(paste("Reading: ", DHS_dataset, sep = ""))
    
    # The DHS dataset of Thurman et al, Nature (2012) contains 8 columns.
    # Column 5-7 coordinates of the distal, non-promoter DHS
    # Column 1-3 contain the connected promoter DHS (= gene).
    DHS_data <- read.delim(DHS_dataset, stringsAsFactors = F, header = F)
    names(DHS_data) <- c("gene_chr", "gene_start", "gene_end", "hgnc_symbol", "dhs_chr","dhs_start","dhs_end", "cor")
    
    DHS_data[,1] <- gsub(pattern = "chr", replacement = "", x = DHS_data[,1])
    DHS_data[,5] <- gsub(pattern = "chr", replacement = "", x = DHS_data[,5])
    
    head(DHS_data)

    # Some distal DHS are located within the gene itself (this may give errors during overlapping)
    DHS_distal <- DHS_data[-which(DHS_data$dhs_start > DHS_data$gene_start-100 & DHS_data$dhs_start < DHS_data$gene_end+100),]
    DHS_distal <- DHS_distal[-which(DHS_distal$dhs_end > DHS_distal$gene_start-100 & DHS_distal$dhs_end < DHS_distal$gene_end+100),]
    
    # Select the genes near the breakpoint junctions
    DHS_affected_genes <- DHS_distal[which(DHS_distal$hgnc_symbol %in% Genes_near_SVs$hgnc_symbol),]
    
    # Count the number of DHS connections per gene
    number_DHS_connections_gene <- data.frame(hgnc_symbol = names(table(DHS_affected_genes$hgnc_symbol)),
                                           Distal_DHS_connections = as.numeric(table(DHS_affected_genes$hgnc_symbol)))
    
    Genes_near_SVs <- merge(Genes_near_SVs, number_DHS_connections_gene, by = "hgnc_symbol", all.x = T)
    
    ## Overlap the "loops" connecting the DHSs and genes with the breakpoints
    DHS_g <- GRanges(seqnames = DHS_affected_genes$dhs_chr, 
                     IRanges(start = ifelse(DHS_affected_genes$dhs_end <= DHS_affected_genes$gene_start,
                                            (DHS_affected_genes$dhs_start+DHS_affected_genes$dhs_end) / 2, 
                                            (DHS_affected_genes$gene_start+DHS_affected_genes$gene_end) / 2),
                             end = ifelse(DHS_affected_genes$dhs_end <= DHS_affected_genes$gene_start,
                                          (DHS_affected_genes$gene_start+DHS_affected_genes$gene_end) / 2, 
                                         (DHS_affected_genes$dhs_start+DHS_affected_genes$dhs_end) / 2)))
    
    Breakpoints_g <- GRanges(seqnames = Input_Breakpoints$Chr, IRanges(start=Input_Breakpoints$Breakpoint, 
                                                                       end = Input_Breakpoints$Breakpoint+1),
                             Patient = Input_Breakpoints$Patient)
  
    olap_DHS_breakpoints <- findOverlaps(DHS_g, Breakpoints_g)
    
    disrupted_DHS_connections <- DHS_affected_genes[queryHits(olap_DHS_breakpoints),]
    disrupted_DHS_connections <- cbind(disrupted_DHS_connections, Input_Breakpoints[subjectHits(olap_DHS_breakpoints),])
    
    # Some genes may be disrupted in multiple different individuals. Therefore filter based on patient id and gene name
    disrupted_DHS_connections$Gene_ID <- paste(disrupted_DHS_connections$Patient, disrupted_DHS_connections$hgnc_symbol, sep = "_")
    
    # Count the number of disrupted DHS connections per gene
    Number_disrupte_DHS <- data.frame(Gene_ID = as.character(names(table(disrupted_DHS_connections$Gene_ID))),
                                      Disrupted_DHS = as.numeric(table(disrupted_DHS_connections$Gene_ID)), stringsAsFactors = F)
    
    Number_disrupte_DHS <- Number_disrupte_DHS[!duplicated(Number_disrupte_DHS$Gene_ID),]
    Genes_near_SVs <- merge(Genes_near_SVs, Number_disrupte_DHS, by = "Gene_ID", all.x = T)
    
    Genes_near_SVs$Distal_DHS_connections[is.na(Genes_near_SVs$Distal_DHS_connections)] <- 0
    Genes_near_SVs$Disrupted_DHS[is.na(Genes_near_SVs$Disrupted_DHS)] <- 0

    # Write the disrupted DHS connections to a new file
    if(output_folder != ""){
      print("# Writing output DHS analysis")
      write.table(x = disrupted_DHS_connections, file = output_file_connections, sep = "\t", quote = F, row.names = F)
      write.table(x = Genes_near_SVs, file = output_file_genes, sep = "\t", quote = F, row.names = F)
    }
  } else {
    print("Re-using previously determined disrupted gene-DHS connections")
    print(paste("Reading: ", output_file_genes))
    Genes_near_SVs <- read.delim(output_file_genes, stringsAsFactors = F)
  }
  
  print("# Finished DHS analysis")
  return(Genes_near_SVs)
}




