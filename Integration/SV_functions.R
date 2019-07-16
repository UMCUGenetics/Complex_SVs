


filter_SVs <- function(input_SVs = SVs,
                       SV_mode = "SVs",
                       limit = distance_limit,
                       Patients_to_exclude = Exclude_patients){
  # The SVs table contains whole affected chromosomes (from start until the end of the chromosome), but we're mainly interested in genes surrounding the SV breakpoint junctions.
  # This function adds fragments to the SV list that are at a maximum distance (limit) from the SV breakpoint junctions.
  

  # Remove the patients in the Exclude_patients group:
  input_SVs <- input_SVs[!(input_SVs$Patient %in% Patients_to_exclude),]
  
  #print(input_SVs)
  output_SVs <- data.frame()
  
  for(patient in unique(input_SVs$Patient)){
    #print(patient)
    SVs_patient <- input_SVs[which(input_SVs$Patient == patient),]

    if(SV_mode == "SVs"){
      for(chromosome in unique(SVs_patient$chr)){
        #print(chromosome)
        
        SVs_patient_chr <- SVs_patient[which(SVs_patient$chr == chromosome),]
        
        SVs_patient_chr$Fragment_start <- SVs_patient_chr$start
        SVs_patient_chr$Fragment_end <- SVs_patient_chr$end
        
        # The first fragment on each affected chromosome runs from basepair one until the first breakpoint, which can be more than the max distance. 
        # Similarly the last fragment goes from the last breakpoint until the end of the chromosome. 
        # These fragments are shortened, so they are at the max distance from the breakpoints.
        if(nrow(SVs_patient_chr) > 1){
          SVs_patient_chr[1,"Fragment_start"] <- SVs_patient_chr[1,"end"] - limit # Select the first breakpoint on the chromosome 
          SVs_patient_chr[nrow(SVs_patient_chr),"Fragment_end"] <- SVs_patient_chr[nrow(SVs_patient_chr),"start"] + limit # Select the last breakpoint on the chromosome 
        }
        output_SVs <- rbind(output_SVs, SVs_patient_chr)
        
      }
      output_SVs$Fragment_start <- ifelse(output_SVs$Fragment_start < 0, 1, output_SVs$Fragment_start) ## Correct for negative values start values 
      
    } else if (SV_mode == "Junctions"){
      
      balanced_data <- SVs_patient[!SVs_patient$ori %in% c("del", "dup"),]
      if(nrow(balanced_data) > 0){
        balanced_svs1 <- data.frame(Patient = patient, 
                                    SV_chr = balanced_data$chr1, 
                                    #SV_start = balanced_data$der1, 
                                    #SV_end = balanced_data$der1, 
                                    SV_ori = balanced_data$ori, 
                                    SV_type = balanced_data$type,
                                    Breakpoint = balanced_data$der1)
        balanced_svs2 <- data.frame(Patient = patient, SV_chr = balanced_data$chr2, 
                                    #SV_start = balanced_data$der2, 
                                    #SV_end = balanced_data$der2, 
                                    SV_ori = balanced_data$ori, 
                                    SV_type = balanced_data$type,
                                    Breakpoint = balanced_data$der2)
        
        
        balanced_svs <- rbind(balanced_svs1, balanced_svs2)
        balanced_svs$SV_start <- balanced_svs$Breakpoint
        balanced_svs$SV_end <- balanced_svs$Breakpoint
      } else {
        balanced_svs <- data.frame(Patient = NULL, SV_chr = NULL,SV_start = NULL, SV_end =NULL, SV_ori = NULL, SV_type = NULL)
      }

      
      cnv_data <- SVs_patient[which(SVs_patient$ori %in% c("del", "dup")),]
      if(nrow(cnv_data) > 0){
        cnvs <- data.frame(Patient = patient, 
                           SV_chr = cnv_data$chr1, 
                           SV_start = cnv_data$der1, 
                           SV_end = cnv_data$der2, 
                           SV_ori = cnv_data$ori,
                           SV_type = cnv_data$type)
      } else {
        cnvs <- data.frame()
      }
      if(nrow(balanced_svs) > 0){
        svs <- rbind(balanced_svs[,c("Patient", "SV_chr", "SV_start", "SV_end", "SV_ori", "SV_type")], cnvs)
      } else {
        svs <- cnvs
      }
      svs$Fragment_start <- svs$SV_start-distance_limit
      
      svs$Fragment_start[which(svs$Fragment_start < 1)] <- 1
      svs$Fragment_end <- svs$SV_end+distance_limit
      
      
      svs$Fragment <- paste(svs$Patient, 1:nrow(svs), sep = "_")
      output_SVs <- rbind(output_SVs, svs)
      
    }
  }
  if(SV_mode == "SVs"){
    output_SVs$Fragment_width <- output_SVs$Fragment_end - output_SVs$Fragment_start
    names(output_SVs)[names(output_SVs) == "width"] <- "SV_width"
    names(output_SVs)[names(output_SVs) == "end"] <- "SV_end"
    names(output_SVs)[names(output_SVs) == "start"] <- "SV_start"
    names(output_SVs)[names(output_SVs) == "chr"] <- "SV_chr"
    names(output_SVs)[names(output_SVs) == "type"] <- "SV_type"
  }
  
  return(output_SVs)
}


Retrieve_breakpoints <- function(input_SVs, SV_mode = "SVs"){
  # This function determines the breakpoint coordinates based on the SV list.
  # A breakpoint is defined as the end of each SV fragment.
  # The last fragment of each chromosome (which goes from the last breakpoint to the end of the chromosome) is skipped (else the end of the chromosome would be classified as a breakpoint).
  output_breakpoints <- data.frame()
  if(SV_mode == "SVs"){
    
    for(patient in unique(input_SVs$Patient)){
      #print(patient)
      SVs_patient <- input_SVs[which(input_SVs$Patient == patient),]
      
        for(chromosome in unique(SVs_patient$SV_chr)){
          #print(chromosome)
          SVs_patient_chr <- SVs_patient[which(SVs_patient$SV_chr == chromosome),]
          
          if(nrow(SVs_patient_chr) > 1){
            breakpoint <- data.frame(Patient = patient, Chr = chromosome, Breakpoint = SVs_patient_chr[,"SV_end"])
            breakpoint <- breakpoint[-nrow(breakpoint),]
          } else {
            breakpoint_start <- data.frame(Patient = patient, Chr = chromosome, Breakpoint = SVs_patient_chr[,"SV_start"])
            breakpoint_end <- data.frame(Patient = patient, Chr = chromosome, Breakpoint = SVs_patient_chr[,"SV_end"])
            breakpoint <- rbind(breakpoint_start,breakpoint_end)
          }
          #print(breakpoint)
          output_breakpoints <- rbind(output_breakpoints, breakpoint)
        }
    }
    } else if (SV_mode == "Junctions"){
      breakpoints1 <- data.frame(Patient = input_SVs$Patient, Chr = input_SVs$SV_chr, Breakpoint = input_SVs$SV_start)
      breakpoints2 <- data.frame(Patient = input_SVs$Patient, Chr = input_SVs$SV_chr, Breakpoint = input_SVs$SV_end)
      breakpoints <- rbind(breakpoints1, breakpoints2)
      breakpoints$BP_ID <- paste(breakpoints$Patient, breakpoints$Chr, breakpoints$Breakpoint, sep = "_")
      output_breakpoints <- breakpoints[!duplicated(breakpoints$BP_ID),]
    }

  return(output_breakpoints)
}



Overlap_Genes_SVs <- function(Patients = levels(SVs$Patient),
                              SVs = all_SVs,
                              Genes = protein_coding_genes,
                              Breakpoints = Breakpoints,
                              Max_Distance_from_BP = 2e6,
                              SV_mode = "SVs"){
  
  print(paste("Determining the genes affected by the SVs (+/- ", Max_Distance_from_BP," bps) in ", length(Patients), " patients", sep = ""))
  
  Genes_SVs_output <- data.frame()
  Genes_G <- GRanges(seqnames = Genes$chromosome_name, IRanges(start = Genes$start_position, end = Genes$end_position))
  
  for(Patient in Patients){
    #print(Patient)
    
    # Select the SVs and Breakpoints of the patient:
    SVs_Patient <- SVs[which(SVs$Patient == Patient),]
    
    if(SV_mode == "SVs"){
      #print(head(SVs_Patient))
      SVs_Patient_G <- GRanges(seqnames = SVs_Patient$SV_chr, IRanges(start = SVs_Patient$Fragment_start, end = SVs_Patient$Fragment_end))
      Breakpoints_Patient <- Breakpoints[which(Breakpoints$Patient == Patient),]
      Breakpoints_Patient_G <- GRanges(seqnames = Breakpoints_Patient$Chr, IRanges(start = Breakpoints_Patient$Breakpoint, end = Breakpoints_Patient$Breakpoint))
      
      #print(SVs_Patient)
      
      # Overlap the genes with the SVs to determine which genes are affected:
      Overlap_SVs_Genes_Patient <- findOverlaps(Genes_G, SVs_Patient_G)
      
      Genes_SV_Patient <- Genes[queryHits(Overlap_SVs_Genes_Patient),]
      Genes_SV_Patient$Patient <- Patient
      Genes_SV_Patient <- cbind(Genes_SV_Patient, SVs_Patient[subjectHits(Overlap_SVs_Genes_Patient), c("Fragment", "Fragment_start", "Fragment_end", "SV_start", "SV_end")])
      
      Genes_SV_Patient$Gene_ID <- paste(as.vector(Genes_SV_Patient$Patient), Genes_SV_Patient$hgnc_symbol, sep = "_")
      
      # Other SV types are added later, therefore change SV_Type to character instead of factor:
      Genes_SV_Patient$SV_type <- as.character(SVs_Patient[subjectHits(Overlap_SVs_Genes_Patient), "SV_type"])
      
      # Truncated genes are overlapped by at least 2 SV fragments and therefore appear multiple times in the list. 
      # The truncated genes will be split according to the fragment they are located on. 
      # The coordinates of these fragmented genes are recalculated to start_position2 and end_position2
      Genes_SV_Patient$start_position2 <- Genes_SV_Patient$start_position
      Genes_SV_Patient$end_position2 <-  Genes_SV_Patient$end_position
      
      # Next the distance from the genes near SVs to the closest breakpoint is determined.
      # The start/end of SV fragments may not overlap the breakpoints when they flank the SVs. 
      # Therefore the distance to breakpoints is calculated and not the distance to the start or end of the SV fragments
      Genes_SV_Patient_G <- GRanges(seqnames = Genes_SV_Patient$chromosome_name, 
                                    IRanges(start = Genes_SV_Patient$start_position2, end = Genes_SV_Patient$end_position2), 
                                    ensembl_gene_id = Genes_SV_Patient$ensembl_gene_id)
      
      # Determine which breakpoint is closest to each gene     
      Distances_Genes_Breakpoints <- distanceToNearest(Genes_SV_Patient_G, Breakpoints_Patient_G)
      Breakpoints_Nearest_Genes <- Breakpoints_Patient_G[subjectHits(Distances_Genes_Breakpoints),]
      
      Genes_SV_Patient2 <- Genes_SV_Patient[queryHits(Distances_Genes_Breakpoints),]
      
      # Add the distance and the nearest breakpoint to the genelist
      Genes_SV_Patient2$Distance_to_BP <- mcols(Distances_Genes_Breakpoints)$distance
      Genes_SV_Patient2$Nearest_BP <- start(Breakpoints_Patient_G[subjectHits(Distances_Genes_Breakpoints),])
      
      # Select the genes that overlap with multiple SV fragments:
      truncated_genes <- Genes_SV_Patient2[which(Genes_SV_Patient2$hgnc_symbol %in% Genes_SV_Patient2$hgnc_symbol[duplicated(Genes_SV_Patient2$ensembl_gene_id)] & Genes_SV_Patient2$Distance_to_BP < 100),]
      
      # Name each fragment of the truncated genes:
      truncated_genes$Gene_ID <- paste(truncated_genes$Gene_ID, truncated_genes$Fragment, sep ="_")
      
      if(nrow(truncated_genes) > 0){
        # First remove the truncated genes from the genelist:
        Genes_SV_Patient2 <- Genes_SV_Patient2[-which(Genes_SV_Patient2$hgnc_symbol %in% truncated_genes$hgnc_symbol),]
        
        # Recalculate  the start and end positions of each gene fragment:
        for(i in 1:nrow(truncated_genes)){
          if(truncated_genes$Fragment_start[i] > truncated_genes$start_position[i] & truncated_genes$Fragment_start[i] < truncated_genes$end_position[i]){
            truncated_genes$start_position2[i] <- truncated_genes$Fragment_start[i]
          }
          if(truncated_genes$Fragment_end[i] > truncated_genes$start_position[i] & truncated_genes$Fragment_end[i] < truncated_genes$end_position[i]){
            truncated_genes$end_position2[i] <- truncated_genes$Fragment_end[i]
          }
        }
        # Add the truncated genes back to the genelist
        truncated_genes$SV_type <- "Truncation"
        Genes_SV_Patient2 <- rbind(Genes_SV_Patient2, truncated_genes)
      }
      
      if(nrow(Genes_SV_Patient2[Genes_SV_Patient2$Distance_to_BP > Max_Distance_from_BP && 
           Genes_SV_Patient2$SV_type != "Deletion" && 
           Genes_SV_Patient2$SV_type != "Duplication",]) > 0){
        Genes_SV_Patient2 <- Genes_SV_Patient2[-which(Genes_SV_Patient2$Distance_to_BP > Max_Distance_from_BP && 
                                                       Genes_SV_Patient2$SV_type != "Deletion" && 
                                                       Genes_SV_Patient2$SV_type != "Duplication"),]
      }
      
      # Some genes may appear more often in the list (usually when they are on insertions/duplications).
      # Remove the duplicated entries. Select based on SV type.
      Genes_SV_Patient2 <- Genes_SV_Patient2[order(Genes_SV_Patient2$SV_type, decreasing = F),] 
      Genes_SV_Patient2 <- Genes_SV_Patient2[!duplicated(Genes_SV_Patient2$Gene_ID),]

      # Add the genes overlapping the SVs in the patient to the complete genelist:
      Genes_SVs_output <- rbind(Genes_SVs_output, Genes_SV_Patient2)
      #print(nrow(Genes_SVs))
    } else if (SV_mode == "Junctions") {
      chr_sizes <- data.frame(chr = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chrX","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr20","chrY","chr19","chr22","chr21","chrM"),
                              size = c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,155270560,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,63025520,59373566,59128983,51304566,48129895,16571))
      
      SVs_g <- GRanges(seqnames = SVs_Patient$SV_chr, IRanges(SVs_Patient$SV_start, SVs_Patient$SV_end), type = SVs_Patient$type)
      
      #all_svs <- c(all_svs,svs_export)
      olap <- distanceToNearest(Genes_G, SVs_g)
      olap_filter <- olap[mcols(olap)$distance < Max_Distance_from_BP,]
      
      overlapping_genes <- cbind(Genes[queryHits(olap_filter),], 
                                 SVs_Patient[subjectHits(olap_filter),], mcols(olap_filter)$distance)
      names(overlapping_genes)[length(names(overlapping_genes))] <- "Distance_to_BP"
      overlapping_genes$Patient <- Patient
      #overview_patient <- summary(overlapping_genes$type)
      overlapping_genes$Gene_ID <- paste(overlapping_genes$Patient, overlapping_genes$hgnc_symbol, sep = "_")
      overlapping_genes$SV_type <- as.character(overlapping_genes$SV_type)
      overlapping_genes$SV_type[which(overlapping_genes$Distance_to_BP < 100 & overlapping_genes$SV_ori != "del" & overlapping_genes$SV_ori != "dup")] <- "Truncation"
      
      overlapping_genes$Fragment_start <- overlapping_genes$SV_start
      overlapping_genes$Fragment_end <- overlapping_genes$SV_end
      
      for(i in 1:nrow(overlapping_genes)){
        if(overlapping_genes$SV_type[i] != "Deletion" & overlapping_genes$SV_type[i] != "Duplication"){
          if(overlapping_genes$transcription_start_site[i] > overlapping_genes$SV_start[i]){
            #print(overlapping_genes[i,])
            end_sv <- min(SVs_Patient$SV_start[SVs_Patient$SV_start > overlapping_genes$SV_start[i] & SVs_Patient$SV_chr == overlapping_genes$SV_chr[i]])
            #print(end_sv)
            if(end_sv == "Inf"){
              overlapping_genes$SV_end[i] <- chr_sizes[gsub(chr_sizes[,1], pattern = "chr", replacement = "") == overlapping_genes$SV_chr[i], 2]
              overlapping_genes$Fragment_end[i] <- overlapping_genes$SV_start[i]+Max_Distance_from_BP
            } else {
              overlapping_genes$SV_end[i] <- end_sv
              overlapping_genes$Fragment_end[i] <- end_sv
              if(overlapping_genes$Fragment_end[i] - overlapping_genes$Fragment_start[i] > Max_Distance_from_BP){
                overlapping_genes$Fragment_end[i] <- overlapping_genes$Fragment_start[i] + Max_Distance_from_BP
              }
            }
              # if geen getal, dan 1
            } else {
              start_sv <- min(SVs_Patient$SV_start[SVs_Patient$SV_start < overlapping_genes$SV_start[i] & SVs_Patient$SV_chr == overlapping_genes$SV_chr[i]])
              #print(start_sv)
              if(start_sv == "Inf"){
                overlapping_genes$SV_start[i] <- 1
                overlapping_genes$Fragment_start[i] <- overlapping_genes$SV_end[i]-Max_Distance_from_BP
              } else {
                overlapping_genes$SV_start[i] <- start_sv
                overlapping_genes$Fragment_start[i] <- start_sv
              }
            }
        }
        overlapping_genes$SV_type[i] <- ifelse(overlapping_genes$Distance_to_BP[i] > 0, "Flanking", overlapping_genes$SV_type[i])
        
        }
      Genes_SVs_output <- rbind(Genes_SVs_output, overlapping_genes)
    }
  }
  
  if(SV_mode == "SVs"){
    # Some (not deleted or duplicated) genes are located further than 2 Mb (Max_Distance_from_BP) from the breakpoints. Remove these genes:
    Genes_SVs_output <- Genes_SVs_output[!(grepl("Normal", Genes_SVs_output$SV_type) & Genes_SVs_output$Distance_to_BP > Max_Distance_from_BP),]
    Genes_SVs_output <- Genes_SVs_output[!(grepl("Inversion", Genes_SVs_output$SV_type) & Genes_SVs_output$Distance_to_BP > Max_Distance_from_BP),]
    #Genes_SVs_output$Patient <- factor(Genes_SVs_output$Patient)
  }
  print(summary(factor(Genes_SVs_output$Patient)))
  
  return(Genes_SVs_output)
}
