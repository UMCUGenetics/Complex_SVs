Run_Phenomatch <- function(phenomatch_folder, # Path to the folder containing the phenomatch program
                           phenomatch_raw_output_folder = paste(phenomatch_folder, "Results/", sep = ""), # The phenomatch program produces raw output files, which are read and processed
                           output_folder, # The final output of this script will be written to this folder (one file per patient and one merged output file containing all patients)
                           output_filename = "Phenomatch_Overview.txt", # Name of the output file containing the merged results for all patients
                           entrez_to_hgnc = TRUE, #Phenomatch uses genenames from HPO, which may not be the same as the ones used in hg19. This uses the hgnc_symbols from genelist, based on entrezgene ids
                           gene_information = "", #only necessary if entrez_to_hgnc = TRUE (needed for the right hgnc_symbols), should contain the data
                           SVs = all_SVs,
                           SV_chr = "SV_chr", # name of the column in the SV dataframe containing the chromosome ids
                           SV_start = "Fragment_start", # name of the column in the SV dataframe containing the start position of the SV
                           SV_end = "Fragment_end", # name of the column in the SV dataframe containing the end position of the SV
                           SV_ID = "Fragment", # name of the column in the SV dataframe containing the IDs of the SVs
                           HPO_obo_file, # path to the hp.obo file
                           HPO_genes_phenotype, # path to the ALL_SOURCES_FREQUENT_FEATURES_genes_to_phenotype file
                           HPO_patients_file, # Path to the file containing the patient IDs and the HPO terms of these patients 
                           overwrite = TRUE,
                           generate_plot = TRUE,
                           patients = ""){
  
  library(reshape2)
  library(ggplot2)
  if(entrez_to_hgnc == TRUE){
    library(biomaRt)
    ensembl_hg19 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # this one is often offline...
  }
  
  print("# Calculating phenomatch for each gene near SV ")
  print(output_folder)
  output_file <- paste(output_folder, output_filename, sep = "")
  dir.create(output_folder, showWarnings = F)
  
  if(file.exists(output_file) == FALSE || overwrite == TRUE ){
    if(file.exists(paste(phenomatch_folder,"bin/phenomatch.jar", sep = ""))){
      print(paste("Reading: ", HPO_genes_phenotype, sep = ""))
      HPO_raw <- read.delim(HPO_genes_phenotype, skip = 1, header = F)
      names(HPO_raw) <- c("entrezgene", "Entrez_gene_name", "HPO_Term", "HPO_Term_ID")
      HPO_raw$entrezgene <- factor(HPO_raw$entrezgene)
      # This file should contain the HPO terms for each patient (should contain a "Patient" and "HPO" (containing HPO IDs) column)
      print(paste("Reading: ", HPO_patients_file, sep = ""))
      HPO_patients <- read.delim(HPO_patients_file)
      
      # The script will only determine the Phenomatch scores of the genes within this distance from the SVs:
      phenomatch_overview_export <- data.frame()
      
      if(patients != ""){
        Patients <- patients
      } else {
        Patients <- as.character(unique(SVs$Patient))
      }
      
      for(patient in Patients){
        output_patient <- paste(output_folder, "Phenomatch_", patient, ".txt", sep = "")
        if(file.exists(output_patient) & overwrite == FALSE){
          print(paste("# Re-using ", output_patient, sep = ""))
          phenomatch_overview <- read.delim(output_patient)
          if(nrow(phenomatch_overview) > 0){
            phenomatch_overview_export_patient <- phenomatch_overview[,c("hgnc_symbol", "phenoMatchScore")]
            phenomatch_overview_export_patient$Patient <- patient
            phenomatch_overview_export_patient$Gene_ID <- paste(patient, phenomatch_overview_export_patient$hgnc_symbol, sep = "_")
            phenomatch_overview_export <- rbind(phenomatch_overview_export, phenomatch_overview_export_patient) 
          }
        } else {
          #print(patient)
          SVs_patient <- SVs[which(SVs$Patient == patient),]
          if(nrow(SVs_patient) > 0){
            HPO_patient <- as.vector(HPO_patients$HPO[which(HPO_patients$Patient == patient)])
            if(length(HPO_patient) > 0){
              HPO_patient <- gsub(",", ";", HPO_patient)
              HPO_patient_split <- unlist(strsplit(HPO_patient, split = ";"))
            }
            # Phenomatch requires a txt file containing the columns: chr - start - end - ID - HPO terms patient
            # One row per SV
            phenomatch_input <- data.frame(chr = paste("chr", SVs_patient[,SV_chr], sep = ""),
                                           start = SVs_patient[,SV_start], 
                                           end = SVs_patient[,SV_end], 
                                           ID = paste(patient, SVs_patient[,SV_ID], sep ="_"),
                                           HPO_Patient = HPO_patient)
            
            # For each gene a Phenomatch score is calculated for all HPO terms of the patient together.
            # In addition, we also want to calculate one score for each phenotypic trait separately 
            # (for example gene A has as score of 19 for HPO:1 and HPO:3 together, and it has a score of 11 for HPO:1 and 12 for HPO:3)
            # Therefore loop over each HPO term seperately:
            for(hpo_term in HPO_patient_split){
              #print(hpo_term)
              
              phenomatch_input_single_HPO <- data.frame(chr = paste("chr", SVs_patient[,SV_chr], sep = ""), 
                                                        start = SVs_patient[,SV_start], 
                                                        end = SVs_patient[,SV_end], 
                                                        ID = paste(patient, SVs_patient[,SV_ID], sep ="_"),
                                                        HPO_Patient = hpo_term)
              
              phenomatch_input <- rbind(phenomatch_input, phenomatch_input_single_HPO)
            }
            
            # The input file for phenomatch will be wrtten here:
            phenomatch_input_file <- paste(phenomatch_folder, "input/", patient, "_phenomatch_input.txt", sep = "")
            
            write.table(file = phenomatch_input_file, x = phenomatch_input, quote = F, sep = "\t", 
                        col.names = F, row.names = F)
            
            phenomatch_output_file <- paste(phenomatch_raw_output_folder, patient, "_phenomatch_output.txt", sep = "")
            
            dir.create(file.path(phenomatch_raw_output_folder), showWarnings = FALSE)
            
            # This is the command to run Phenomatch.
            phenomatch_command <- paste("java -jar ", 
                                        phenomatch_folder,"bin/phenomatch.jar -i ", 
                                        phenomatch_input_file, " -g ", 
                                        phenomatch_folder, "data/knownGene.txt.entrez_id.tab.unique -O ", 
                                        HPO_obo_file, " -a ", 
                                        HPO_genes_phenotype, " -o ", 
                                        phenomatch_output_file, sep = "")
            print(phenomatch_command)
            system(phenomatch_command)
            
            # Phenomatch automatically adds ".overlapped_genes.txt" to the output filename
            phenomatch_file <- paste(phenomatch_output_file, ".overlapped_genes.txt", sep = "")
            
            # Read the Phenomatch output
            phenomatch_patient <- read.delim(phenomatch_file, stringsAsFactors = F)
            
            # Obtain the phenomatch score for the full phenotype (all HPO terms together) by removing the single HPO terms; except when here is only one HPO-term 
            if(length(HPO_patient_split) > 1){
              phenomatch_patient_full <- phenomatch_patient[!as.vector(phenomatch_patient$phenotypes) %in% HPO_patient_split,]
            }else{
              phenomatch_patient_full <- phenomatch_patient
            }
            
            phenomatch_overview <- data.frame(gene_symbol = phenomatch_patient_full$gene_symbol, 
                                              #entrezgene = phenomatch_patient_full$entrezgene,
                                              phenoMatchScore = phenomatch_patient_full$phenoMatchScore, stringsAsFactors = F)
            phenomatch_overview <- phenomatch_overview[!duplicated(phenomatch_overview$gene_symbol),]
            phenomatch_overview <- phenomatch_overview[which(phenomatch_overview$phenoMatchScore > 0),]
            
            # Add the phenomatch scores for single HPO terms as columns to the phenomatch_overview
            for(hpo_term in HPO_patient_split){
              #print(hpo_term)
              
              phenomatch_single_HPO <- phenomatch_patient[phenomatch_patient$phenotypes == hpo_term,]
              
              phenomatch_single_HPO_filter <- phenomatch_single_HPO[,c("gene_symbol","phenoMatchScore")]
              phenomatch_single_HPO_filter <- phenomatch_single_HPO_filter[!duplicated(phenomatch_single_HPO_filter$gene_symbol),]
              phenomatch_single_HPO_filter <- phenomatch_single_HPO_filter[which(phenomatch_single_HPO_filter$phenoMatchScore > 0),]
              names(phenomatch_single_HPO_filter)[2] <- hpo_term
              phenomatch_overview <- merge(phenomatch_overview, phenomatch_single_HPO_filter, all.x = T)
            }
            
            # Translate the HPO ids to HPO terms:
            all_hpos <- HPO_raw[,c("HPO_Term", "HPO_Term_ID")]
            all_hpos <- all_hpos[!duplicated(all_hpos$HPO_Term_ID),]
            
            hpos <- data.frame(HPO_Term_ID = HPO_patient_split)
            hpos_patient <- merge(hpos, all_hpos, by = "HPO_Term_ID")
            
            hpo <- c(names(phenomatch_overview))
            hpo[1] <- "HPO_ID"
            phenomatch_overview$gene_symbol <- as.vector(phenomatch_overview$gene_symbol)
            phenomatch_overview <- rbind(phenomatch_overview, hpo)
            
            for(names in names(phenomatch_overview)){
              #print(names)
              
              hpo_id <- hpos_patient[hpos_patient$HPO_Term_ID == names,]
              if(nrow(hpo_id) > 0){
                #print(hpo_id$HPO_Term)
                names(phenomatch_overview)[which(names(phenomatch_overview) == names)] <- as.vector(hpo_id$HPO_Term)
                
              }
            }
            #print(phenomatch_overview[nrow(phenomatch_overview),])
            
            phenomatch_overview <- phenomatch_overview[which(phenomatch_overview$gene_symbol != "HPO_ID"),]
            phenomatch_overview$phenoMatchScore <- as.numeric(phenomatch_overview$phenoMatchScore)
            phenomatch_overview$gene_symbol <- factor(phenomatch_overview$gene_symbol, levels = phenomatch_overview$gene_symbol[order(phenomatch_overview$phenoMatchScore)])
            
            # Phenomatch uses the genenames from HPO, which are in some cases not the same as the (old) hg19 HGNC symbols we use.
            # Therefore translate the phenomatch genenames to hg19 HGNC symbols:
            if(entrez_to_hgnc == TRUE){
              entrez_ids <- HPO_raw[!duplicated(HPO_raw$entrezgene), c(1,2)]
              names(entrez_ids) <- c("entrezgene", "gene_symbol")
              # add the entrezgene ids to phenomatch_overview
              phenomatch_overview <- merge(phenomatch_overview, entrez_ids, by = "gene_symbol", all.x = T)
              # add the hgnc_symbols from the genelist to the phenomatch_overview based on entrezgene ids
              phenomatch_overview <- merge(phenomatch_overview, gene_information[,c("entrezgene", "hgnc_symbol")], by = "entrezgene", all.x = T)
            }
            write.table(file = paste(output_folder, "Phenomatch_", patient, ".txt", sep = ""), x = phenomatch_overview[,c("hgnc_symbol", "phenoMatchScore", names(phenomatch_overview)[4:(ncol(phenomatch_overview)-1)])], row.names = F, sep = "\t", quote = F)
          }
          
          if(nrow(phenomatch_overview) > 0){
            if(entrez_to_hgnc == TRUE){
              phenomatch_overview_export_patient <- phenomatch_overview[,c("hgnc_symbol", "phenoMatchScore")]
            } else {
              phenomatch_overview_export_patient <- phenomatch_overview[,c("gene_symbol", "phenoMatchScore")]
              names(phenomatch_overview_export_patient)[names(phenomatch_overview_export_patient) == "gene_symbol"] <- "hgnc_symbol"
            }
            
            phenomatch_overview_export_patient$Patient <- patient
            
            phenomatch_overview_export_patient$Gene_ID <- paste(patient, phenomatch_overview_export_patient$hgnc_symbol, sep = "_")
            
            phenomatch_overview_export <- rbind(phenomatch_overview_export, phenomatch_overview_export_patient) 
            
            #print(head(phenomatch_overview))
            
            if(generate_plot == TRUE){
              phenomatch_overview_plot <- phenomatch_overview[,-which(names(phenomatch_overview) %in% c("entrezgene", "hgnc_symbol"))]
              phenomatch_overview_melted <- melt(phenomatch_overview_plot, id.vars = "gene_symbol")
              phenomatch_overview_melted$value <- as.numeric(phenomatch_overview_melted$value)
              phenomatch_overview_melted$Score <- phenomatch_overview_melted$value
              phenomatch_overview_melted$Score <- ifelse(phenomatch_overview_melted$Score > 30, 30, phenomatch_overview_melted$Score)
              
              phenomatch_overview_melted$value <- round(phenomatch_overview_melted$value, 1)
              
              names(phenomatch_overview_melted) <- c("Gene", "Phenotype", "phenoMatchScore", "Score")
              
              output_height <- nrow(phenomatch_overview) / 3 + 2.5
              
              #print(output_height)
              output_width <- ncol(phenomatch_overview) / 2 + 2
              #print(output_width)
              
              data_table <- ggplot(phenomatch_overview_melted, aes(x = Phenotype, y = Gene, fill = Score), label = phenoMatchScore) + geom_tile(lwd = 0.8, colour  = "black") +
                geom_text(size = 3.5, label = phenomatch_overview_melted$phenoMatchScore) + 
                scale_fill_gradient(low ="white", high = "red", limits=c(0,30)) +
                scale_x_discrete(position = "top") + 
                theme_minimal() + 
                theme(axis.text.x = element_text(angle = 60, hjust = 0), axis.text.y = element_text(face = "italic"), plot.title = element_text(hjust = 0.5)) + ggtitle(label = patient)
              
              #print(head(phenomatch_overview_melted))
              
              pdf(file = paste(output_folder, "Phenomatch_", patient, ".pdf", sep = ""), height = output_height, width = output_width)
              
              print(data_table)
              dev.off()
              
              print(paste("Exported plot to: ",  output_folder, "Phenomatch_", patient, ".pdf", sep = ""))
            }
          }
        }
      }
      print(paste("Writing: ", output_folder, "Phenomatch_Overview.txt", sep = ""))
      write.table(file = paste(output_folder, "Phenomatch_Overview.txt", sep = ""), x = phenomatch_overview_export, row.names = F, sep = "\t", quote = F)
      
      return(phenomatch_overview_export)
    } else {
      print(paste("! Error: ", phenomatch_folder,"bin/phenomatch.jar NOT found!", sep = ""))
      stop()
    }
  } else {
    print(paste("Reading: ", output_folder, "Phenomatch_Overview.txt", sep = ""))
    
    phenomatch_overview_export <- read.delim(paste(output_folder, "Phenomatch_Overview.txt", sep = ""), stringsAsFactors = F)
    return(phenomatch_overview_export)
  }
} 


# This function counts the number of overlaps between the HPO terms of a patient with the HPO terms associated with each gene. 
Link_HPO_Genes_to_Patients <- function(Genes = Genes_SVs,
                                       Patients = unique(Genes_SVs$Patient),
                                       HPO_Terms_Patients = HPO_Patients_file,
                                       HPO_file = HPO_genes_phenotype,
                                       Phenomatch_output_folder,  
                                       Phenomatch_threshold = 5){
  
  print(paste("# Calculating the number of phenomatch hits per gene per patient (PhenomatchScore > ", Phenomatch_threshold, ")", sep = ""))
  
  HPO_Patients <- read.delim(HPO_Terms_Patients, stringsAsFactors = F)
  HPO_Terms <- read.delim(HPO_file, skip = 1, header = F)
  names(HPO_Terms) <- c("entrezgene", "Entrez_gene_name", "HPO_Term", "HPO_Term_ID")
  HPO_Terms$entrezgene <- factor(HPO_Terms$entrezgene)
  Genes_SVs_HPO <- data.frame()
  
  for(Patient in Patients){
    #print(Patient)
    
    # Select the genes near SVs for the patient
    Genes_Patient <- Genes[which(Genes$Patient == Patient),]
    Genes_Patient <- Genes_Patient[which(Genes_Patient$entrezgene != ""),]
    Genes_Patient$entrezgene <- factor(Genes_Patient$entrezgene)
    
    # Select the HPO terms of the patients
    HPO_Patient <- HPO_Patients[which(HPO_Patients$Patient == Patient),]
    HPO_Patient$HPO <- as.vector(HPO_Patient$HPO)
    
    HPO_Patient2 <- unlist(strsplit(HPO_Patient$HPO,  split = ","))
    number_HPO_Patient <- length(HPO_Patient2)
    #print(number_HPO_Patient)
    
    if(number_HPO_Patient > 0){
      phenomatch_file <- paste(Phenomatch_output_folder, "Phenomatch_", Patient, ".txt", sep = "")
      phenomatch_patient <- read.delim(phenomatch_file, check.names = F)
      
      if(nrow(phenomatch_patient) > 0){
        phenomatch <- as.matrix(phenomatch_patient[,-c(1,2)])
        phenomatch[phenomatch < Phenomatch_threshold] <- 0
        phenomatch[phenomatch >= Phenomatch_threshold] <- 1
        phenomatch[is.na(phenomatch)] <- 0
        Phenomatch_hits_patient <- data.frame(Patient = Patient, hgnc_symbol = phenomatch_patient$hgnc_symbol, Phenomatches = rowSums(phenomatch))
      } else {
        Phenomatch_hits_patient <- data.frame(Patient = Patient, hgnc_symbol = NA, Phenomatches = 0)
      }
      Genes_SVs_HPO_patient <- merge(Genes_Patient[,c("hgnc_symbol", "Patient")], Phenomatch_hits_patient, by = c("hgnc_symbol", "Patient"), all.x = T)
      Genes_SVs_HPO_patient$Phenomatches[is.na(Genes_SVs_HPO_patient$Phenomatches)] <- 0
      Genes_SVs_HPO_patient$HPO_Patient <- number_HPO_Patient
      Genes_SVs_HPO <- rbind(Genes_SVs_HPO, Genes_SVs_HPO_patient)
    } else {
      print(paste("Patient ", Patient, " does not have any HPO terms", sep = ""))
    } 
  }
  #print("# Done ")
  return(Genes_SVs_HPO)
}


# This function selects the proposed phenodrivers of each patient and determines how many of the HPO terms of these patients are associated with these phenodrivers.
patient_summary <- function(genes_svs,
                            phenomatch_output_folder = "",
                            phenomatch_threshold = 4,
                            conclusions = c("Partially", "Largely"),
                            conclusion_thresholds = c(0.2, 0.6)){
  
  print(paste("# Generating driver overview for ", length(unique(genes_svs$Patient)), " patients.", sep = ""))
  driver_match_overview <- data.frame()
  for(patient in unique(genes_svs$Patient)){
    #print(patient)
    genes_svs_patient <- genes_svs[genes_svs$Patient == patient,]
    phenomatch_file <- paste(phenomatch_output_folder, "Phenomatch_", patient, ".txt", sep = "")
    phenomatch_patient <- read.delim(phenomatch_file, stringsAsFactors = F)
    
    # only select the genes that are affected by the SVs ("phenodrivers"):
    phenomatch_drivers <- phenomatch_patient[phenomatch_patient$hgnc_symbol %in% genes_svs_patient$hgnc_symbol[which(genes_svs_patient$Classification != "No_driver" & 
                                                                                                                       genes_svs_patient$Classification != "T3")],]
    
    # count the number of HPO terms (= number of columns - first 2 columns)
    hpo_terms_patient <- (ncol(phenomatch_drivers)-2)
    
    total_phenomatch_patient <- sum(phenomatch_drivers$phenoMatchScore)
    
    phenomatch_drivers2 <- as.matrix(phenomatch_drivers[,-c(1:2)])
    phenomatch_drivers2[phenomatch_drivers2 < phenomatch_threshold] <- 0
    phenomatch_drivers2[phenomatch_drivers2 > phenomatch_threshold] <- 1
    #colSums(phenomatch_drivers2)
    
    # Determine how many of the patient's HPO terms have at least one gene associated with it (= rowSums > 0):
    if(ncol(phenomatch_drivers2) > 3){
      HPO_matches <- length(which(colSums(phenomatch_drivers2) > 0)) / (ncol(phenomatch_drivers2))
    } else {
      HPO_matches <- length(which(phenomatch_drivers2 > 0)) / (ncol(phenomatch_drivers2))
      HPO_matches <- ifelse(HPO_matches > 1, 1, HPO_matches)
    }
    
    # Add a letter to the gene name indicating if a gene is deleted, duplicated, truncated or flanked by the SVs
    genes_svs_patient$Effect <- ifelse(genes_svs_patient$SV_type %in% c("Deletion", "Duplication","Insertion","Truncation"), "Direct", "Indirect")
    
    genes_svs_patient$gene_symbol <- ifelse(genes_svs_patient$SV_type == "Truncation", paste(genes_svs_patient$hgnc_symbol, " (T)", sep = ""), 
                                            paste(genes_svs_patient$hgnc_symbol, " (F)", sep = ""))
    genes_svs_patient$gene_symbol <- ifelse(genes_svs_patient$SV_type == "Duplication" | genes_svs_patient$SV_type == "Insertion", paste(genes_svs_patient$hgnc_symbol, " (D)", sep = ""), genes_svs_patient$gene_symbol)
    genes_svs_patient$gene_symbol <- ifelse(genes_svs_patient$SV_type == "Deletion", paste(genes_svs_patient$hgnc_symbol, " (d)", sep = ""), genes_svs_patient$gene_symbol)
    
    # collect the candidate drivers per classification tier
    T1_direct <- genes_svs_patient[which(genes_svs_patient$Classification == "T1" & genes_svs_patient$Effect == "Direct"),]
    T1_indirect <- genes_svs_patient[which(genes_svs_patient$Classification == "T1" & genes_svs_patient$Effect == "Indirect"),]
    T2_direct <- genes_svs_patient[which(genes_svs_patient$Classification == "T2" & genes_svs_patient$Effect == "Direct"),]
    T2_indirect <- genes_svs_patient[which(genes_svs_patient$Classification == "T2" & genes_svs_patient$Effect == "Indirect"),]
    T3_direct <- genes_svs_patient[which(genes_svs_patient$Classification == "T3" & genes_svs_patient$Effect == "Direct"),]
    T3_indirect <- genes_svs_patient[which(genes_svs_patient$Classification == "T3" & genes_svs_patient$Effect == "Indirect"),]
    
    patient_conclusion <- "VUS"
    for(i in 1:length(conclusions)){
      patient_conclusion <- ifelse(HPO_matches >= conclusion_thresholds[i], conclusions[i], patient_conclusion)
    }
    
    driver_match_patient <- data.frame(Patient = patient, 
                                       HPO_match = round(HPO_matches,2), 
                                       HPO_terms = hpo_terms_patient, 
                                       phenomatch  = total_phenomatch_patient,
                                       T1_direct = paste(unique(T1_direct$gene_symbol[order(T1_direct$Phenomatches_high, decreasing = T)]), collapse = ","),
                                       T1_indirect = paste(unique(T1_indirect$gene_symbol[order(T1_indirect$Phenomatches_high, decreasing = T)]), collapse = ","),
                                       T2_direct = paste(unique(T2_direct$gene_symbol[order(T2_direct$Phenomatches_high, decreasing = T)]), collapse = ","),
                                       T2_indirect = paste(unique(T2_indirect$gene_symbol[order(T2_indirect$Phenomatches_high, decreasing = T)]), collapse = ","),
                                       T3_direct = paste(unique(T3_direct$gene_symbol[order(T3_direct$Phenomatches_high, decreasing = T)]), collapse = ","),
                                       T3_indirect = paste(unique(T3_indirect$gene_symbol[order(T3_indirect$Phenomatches_high, decreasing = T)]), collapse = ","),
                                       Conclusion = patient_conclusion)
    driver_match_overview <- rbind(driver_match_overview, driver_match_patient)
  }
  return(driver_match_overview)
}
