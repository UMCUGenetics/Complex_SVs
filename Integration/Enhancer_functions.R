## First determine in which cell types each gene is highest expressed
## Subsequently select the TADs and the enhancers of these cell types
## Finally determine the loss of enhancers

# For each gene the 5 cell types with the highest expression of the gene are selected 
# RNA_expression_data should contain ensembl gene ids as row.names and celltypes as col.names
select_highest_expression_gene <- function(RNA_expression_data, 
                                           number_celltypes = 5,
                                           metadata = encode_metadata,
                                           ensembl_gene_ids, # these are the ids for the genes that need expression info
                                           column_ensembl_ids = 1, # This column should contain the ensembl gene ids and will be used as row.names
                                           first_column_expression_data = 4, # all columns before this column will be removed
                                           output = "data.frame"){
  
  #print(paste("# Reading :", RNA_expression_file, sep = ""))
  #RNA_expression_data <- read.delim(RNA_expression_file)
  row.names(RNA_expression_data) <- RNA_expression_data[,column_ensembl_ids]
  RNA_expression_data <- RNA_expression_data[,-c(1:(first_column_expression_data)-1)]
  if(output == "list"){
    celltypes_gene_expression <- list()
    for(gene in ensembl_gene_ids){
      rna_gene <- RNA_expression_data[gene,]
      celltypes_gene_expression[[gene]] <- data.frame("V1" =names(rna_gene[,order(rna_gene, decreasing = T)])[1:number_celltypes])
      #celltypes_gene_expression[[gene]] <- merge(celltypes_gene_expression[[gene]], metadata[,c("V1", "V4")], by = "V1", all.x = T)
    }
  } else if (output ==  "data.frame"){
    celltypes_gene_expression <- data.frame()
    
    for(gene in ensembl_gene_ids){
      rna_gene <- RNA_expression_data[gene,]
      
      celltypes_gene <- data.frame(ensembl_gene_id = gene, "Cell_types" =names(rna_gene[,order(rna_gene, decreasing = T)])[1:number_celltypes])
      celltypes_gene_expression <- rbind(celltypes_gene_expression, celltypes_gene)
    }
  }
  return(celltypes_gene_expression)
}


determine_enhancer_loss_gene <- function(Genes_SVs,
                                         RNA_expression_data,
                                         input_TADs,
                                         folder_18_state_enhancers,
                                         folder_15_state_enhancers){
  
  enhancer_loss_overview <- data.frame()
  
  # First determine in which cell types each gene is highest expressed
  print("# Select the cell types with highest gene expression for each gene")
  celltypes_gene_expression <- select_highest_expression_gene(ensembl_gene_ids = unique(Genes_SVs$ensembl_gene_id),
                                                              RNA_expression_data = RNA_expression_data,
                                                              number_celltypes = 3, output = "data.frame")
  
  celltypes_gene_expression2 <- merge(Genes_SVs[,c("Gene_ID","ensembl_gene_id","chromosome_name","transcription_start_site", "SV_start", "SV_end")],
                                      celltypes_gene_expression, by = "ensembl_gene_id", all.x = T)
  
  # Select the TADs of cell types with gene expression for each gene
  print(paste("# Reading ", input_TADs, sep = ""))
  TADs <- read.delim(input_TADs, stringsAsFactors = F)
  TADs$TAD_chr <- gsub("chr", "", TADs$TAD_chr)
  TADs_per_gene <- data.frame()
  
  enhancer_list <- list()
  
  print("# Determining enhancer loss per cell type")
  for(cell_type_HiC in unique(celltypes_gene_expression2$Cell_types)){
    print(cell_type_HiC)
    
    # Select the TAD for each gene 
    gene_expression_cell <- celltypes_gene_expression2[celltypes_gene_expression2$Cell_types == cell_type_HiC,]
    gene_expression_cell_g <- GRanges(seqnames = gene_expression_cell$chromosome_name, IRanges(start = gene_expression_cell$transcription_start_site, end = gene_expression_cell$transcription_start_site))
    TADs_celltype <- TADs[which(TADs$TAD_cell == cell_type_HiC),]
    TADs_cell_type <- GRanges(seqnames = TADs_celltype$TAD_chr, IRanges(start = TADs_celltype$TAD_start, end = TADs_celltype$TAD_end))
    olap_genes_TADs <- findOverlaps(gene_expression_cell_g, TADs_cell_type)
    
    # Add the TAD coordinates to the affected genes
    Gene_TADs <- cbind(gene_expression_cell[queryHits(olap_genes_TADs),], TADs_celltype[subjectHits(olap_genes_TADs),c("TAD_start", "TAD_end")])
    
    # Determine the portion of the disrupted TAD containing the gene (=remaining TAD)
    Gene_TADs$Remaining_TAD_start <- ifelse(Gene_TADs$SV_start < Gene_TADs$TAD_start, Gene_TADs$TAD_start, Gene_TADs$SV_start)
    Gene_TADs$Remaining_TAD_end <- ifelse(Gene_TADs$SV_end < Gene_TADs$TAD_end, Gene_TADs$SV_end, Gene_TADs$TAD_end)
    
    # In some cases the TAD is not affected by the SV, and therefore the TAD and SV do not overlap. In these cases there is no loss of enhancers. 
    # Set the coordinates of the "Remaining_TAD" to the coordinates of the entire TAD
    Gene_TADs$Remaining_TAD_start[which(Gene_TADs$Remaining_TAD_end < Gene_TADs$Remaining_TAD_start)] <- Gene_TADs$TAD_start[which(Gene_TADs$Remaining_TAD_end < Gene_TADs$Remaining_TAD_start)]
    Gene_TADs$Remaining_TAD_end[which(Gene_TADs$Remaining_TAD_end < Gene_TADs$Remaining_TAD_start)] <- Gene_TADs$TAD_end[which(Gene_TADs$Remaining_TAD_end < Gene_TADs$Remaining_TAD_start)]
    
    # Translate the cell types from the HiC/TAD data to the cell type ID used by Roadmap/Encode
    celltype_roadmap  <- celltype_metadata[celltype_metadata$HiC_ID == cell_type_HiC, "Roadmap_ID"]
    
    # Read the enhancer data for the cell type
    if(celltype_roadmap %in% names(enhancer_list) == FALSE){
      # read the 18-state enhancers if it exists, else read 15-state enhancers
      if(length(list.files(folder_18_state_enhancers)[grep(celltype_roadmap, x = list.files(folder_18_state_enhancers))] > 0)){
        enhancer_file <- paste(folder_18_state_enhancers, list.files(folder_18_state_enhancers)[grep(celltype_roadmap, x = list.files(folder_18_state_enhancers))], sep = "")
        print(paste("Reading: ", enhancer_file, sep = ""))
        enhancers_cell_type <- read.delim(enhancer_file, stringsAsFactors = F)
        enhancers_g <- GRanges(seqnames = gsub(enhancers_cell_type[,1], pattern = "chr", replacement = ""), IRanges(start = enhancers_cell_type[,2], end = enhancers_cell_type[,3]))
      } else {
        enhancer_file <- paste(folder_15_state_enhancers, list.files(folder_15_state_enhancers)[grep(celltype_roadmap, x = list.files(folder_15_state_enhancers))], sep = "")
        print(paste("Reading: ", enhancer_file, sep = ""))
        enhancers_cell_type <- read.delim(enhancer_file, stringsAsFactors = F)
        enhancers_g <- GRanges(seqnames = gsub(enhancers_cell_type[,1], pattern = "chr", replacement = ""), IRanges(start = enhancers_cell_type[,2], end = enhancers_cell_type[,3]))
      }
      # add the enhancers of the celltype to the list so they wont have to be re-read later
      enhancer_list[[celltype_roadmap]] <- enhancers_g
    } else {
      enhancers_g <- enhancer_list[[grep(pattern = celltype_roadmap, x = names(enhancer_list))]]
    }
    
    # Count the number of ernhancers per TAD
    Gene_TADs$Enhancers_TAD <- countOverlaps(GRanges(seqnames = Gene_TADs$chromosome_name, IRanges(start = Gene_TADs$TAD_start, end = Gene_TADs$TAD_end)), enhancers_g)
    
    # Count the enhancers that are located within the part of the disrupted TAD containing the gene 
    Gene_TADs$Enhancers_remaining <- countOverlaps(GRanges(seqnames = Gene_TADs$chromosome_name, 
                                                        IRanges(start = Gene_TADs$Remaining_TAD_start, end = Gene_TADs$Remaining_TAD_end)), enhancers_g)
    Gene_TADs$Enhancers_lost <- (Gene_TADs$Enhancers_TAD - Gene_TADs$Enhancers_remaining) / Gene_TADs$Enhancers_TAD * 100
    enhancer_loss_overview <- rbind(enhancer_loss_overview, Gene_TADs) 
  }
  return(enhancer_loss_overview)
}
