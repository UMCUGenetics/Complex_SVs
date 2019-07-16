
# This function adds RNA-seq RPKM values and differential expression results to each gene in the Genes_SVs table
# The files should contain ensembl_gene_ids as row.names and a column for each individual
Add_RNA_RPKM_DE <- function(RNA_RPKM_files,
                            RNA_DE_Folder,
                            Genes_SVs){
  
  output <- data.frame()
  
  for(RPKM_file in RNA_RPKM_files){
    #print(RPKM_file)
    RNA_RPKMs_raw <- read.table(RPKM_file, check.names = F, header = T)
    
    # The RNA-seq files have names different from the patient IDs. Translate the RNA-seq filenames to patient ids:
    RNA_overview <- data.frame(RNA_files = names(RNA_RPKMs_raw), stringsAsFactors = F)
    # Remove the -RNA-1 suffices
    RNA_overview$RNA_files <- gsub("-RNA-[[:digit:]]$", "", RNA_overview$RNA_files)
    # remove family IDs 
    RNA_overview$Patient_IDs <- gsub("MP[0-9][0-9]-", "MP", RNA_overview$RNA_files)
    # The RNA-seq data files of patient UTR22 are called "Pregno-child-blood", change name to UTR22
    RNA_overview$Patient_IDs[RNA_overview$RNA_files == "Pregno-child-blood"] <- "UTR22"

    for(Patient in unique(RNA_overview$Patient_IDs)){
      #print(Patient)

      if(Patient %in% unique(Genes_SVs$Patient)){
        
        # Select the DESeq2 output file for the patient
        RNA_DE_file_patient <- list.files(RNA_DE_Folder, full.names = T)[grep(RNA_overview$RNA_files[RNA_overview$Patient_IDs == Patient], 
                                                                              list.files(RNA_DE_Folder))]
             
        if(length(RNA_DE_file_patient) == 1){
          print(paste("Reading ", RNA_DE_file_patient, sep = ""))
          RNA_DE_patient <- read.table(RNA_DE_file_patient)
          names(RNA_DE_patient) <- paste("RNA_", names(RNA_DE_patient), sep = "")
          RNA_DE_patient$ensembl_gene_id <- row.names(RNA_DE_patient)

        } else {
          print(paste("! Error loading RNA DE results patient ", Patient, " !", sep = ""))
        }
        
        # Add the RPKMs for the patient and means of the other individuals (Controls) to the Genes_SVs of the patient
        RNA_RPKMs_patient <- data.frame(ensembl_gene_id = row.names(RNA_RPKMs_raw),
                                        RNA_RPKM = rowMeans(RNA_RPKMs_raw[grep(pattern = RNA_overview$RNA_files[RNA_overview$Patient_IDs == Patient], names(RNA_RPKMs_raw))]),
                                        RNA_Controls_RPKM = rowMeans(RNA_RPKMs_raw[-grep(pattern = RNA_overview$RNA_files[RNA_overview$Patient_IDs == Patient], names(RNA_RPKMs_raw))]))
          
        Genes_SVs_patient <- Genes_SVs[which(Genes_SVs$Patient == Patient),]
                
        if(nrow(Genes_SVs_patient) > 0){
          Genes_SVs_patient <- merge(Genes_SVs_patient, RNA_RPKMs_patient, by = "ensembl_gene_id", all.x = T)
                    
          if(nrow(RNA_DE_patient) > 0){
            Genes_SVs_patient <- merge(Genes_SVs_patient, RNA_DE_patient[,c("ensembl_gene_id","RNA_log2FoldChange", "RNA_pvalue", "RNA_padj")], 
                                       by = "ensembl_gene_id", all.x = T)
          }
          
          if(nrow(output) > 0){
            output <- rbind(output, Genes_SVs_patient) 
          } else {
            output <- Genes_SVs_patient
          }
        } else {
         print(paste("# Genes SV data not found for patient: ", Patient, sep = ""))
        }
      } else {
         print(paste("# ",Patient, " not in Genes SVs", sep = ""))
      }
    }
  }
  return(output)
}

## Function to add truncated fragment RNA data; Data replaces Log2FoldChange, RNA_padj and RNA_pvalue and RNA_RPKMs. 
Add_Trun_RNA <- function(RNA_Truncated_DE_Folders, Subset_Folder, RNA_Output){
  
  for(Truncated_DE_folder in RNA_Truncated_DE_Folders){
    #Select here for which patients the truncated genes should be added; this depends on which subset it is calculated
    Subset <- gsub("/DE", "", gsub(".*DE_", "", Truncated_DE_folder))
    Subset <- paste(Subset_Folder, "Truncated_Genes_", Subset, ".txt", sep = "")
    SubsetGenes <- read.delim(Subset, header = T, stringsAsFactors = F)
    SubsetGenessplit <- strsplit(as.character(SubsetGenes$Gene_ID), "_", fixed = T)
    SubsetGenessplit <- data.frame(do.call(rbind, SubsetGenessplit), SubsetGenes$ensembl_gene_id)
    SubsetGenessplit$ID <- paste(SubsetGenessplit$SubsetGenes.ensembl_gene_id, "_", SubsetGenessplit$X3, sep = "")
    
    ## LOAD DATA 
    DE_Files <- list.files(Truncated_DE_folder, full.names = F)[grep(".bed", list.files(Truncated_DE_folder, full.names = F))]
    
    #RNA_patients <- gsub(".*_", "", gsub("_DE.*", "", DE_Files))
    #print(RNA_patients)
    
    # The RNA-seq files have names different from the patient IDs. Translate the RNA-seq filenames to patient ids:
    RNA_patients <- data.frame(RNA_files = gsub(".*_", "", gsub("_DE.*", "", DE_Files)), stringsAsFactors = F)
    # remove family IDs 
    RNA_patients$Patient_IDs <- gsub("MP[0-9][0-9]-", "MP", RNA_patients$RNA_files)
    # The RNA-seq data files of patient UTR22 are called "Pregno-child-blood", change name to UTR22
    RNA_patients$Patient_IDs[RNA_patients$RNA_files == "Pregno-child-blood"] <- "UTR22"
    
    #print(RNA_patients)
    
    # First add the differential expression data of the truncated genes to the RNA_Output
    for(patient in unique(RNA_patients$RNA_files)){
      
      Bedfile <- DE_Files[grep(patient, DE_Files)]
      if(length(Bedfile) == 1){
        
        Trun_bed <- read.delim(paste(Truncated_DE_folder, "/", Bedfile, sep = ""), check.names = F)
        names(Trun_bed) <- paste("RNA_", names(Trun_bed), sep = "")
        
        #Select here only for the truncated genes
        Truncated_genes <- Trun_bed[grep("ENSG[0-9]{11}_", row.names(Trun_bed)),]
        Truncated_genes$ensembl_gene_id <- row.names(Truncated_genes)
        
        TruncatedPatient <- SubsetGenessplit[SubsetGenessplit$X1 == patient,]
        
        if(patient == "Pregno-child-blood"){
          TruncatedPatient <- SubsetGenessplit[SubsetGenessplit$X1 == "Pregno",]
        }
        
        if(nrow(TruncatedPatient) == 0){
          #print(paste("No truncated genes added for", patient, sep = " "))
          next
        }
        
        TruncatedData <- Truncated_genes[Truncated_genes$ensembl_gene_id %in% TruncatedPatient$ID,]
        TruncatedData <- merge(TruncatedData, TruncatedPatient, by.x = "ensembl_gene_id", by.y = "ID")
        TruncatedData$Gene_ID <- paste(RNA_patients$Patient_IDs[RNA_patients$RNA_files == TruncatedData$X1], 
                                       "_", TruncatedData$X2, "_", TruncatedData$X3, sep = "")
        
        if(patient == "Pregno-child-blood"){
          TruncatedData$Gene_ID <- paste(RNA_patients$Patient_IDs[RNA_patients$RNA_files == "Pregno-child-blood"], 
                                         "_", TruncatedData$X2, "_", TruncatedData$X3, sep = "")
        }                               
                                       
        print(paste("Truncated gene expression data added for ", patient , sep = " "))
        
        for(i in 1:nrow(TruncatedData)){
          for(k in 1:nrow(RNA_Output)){
            if(TruncatedData$Gene_ID[i] == RNA_Output$Gene_ID[k]){
              RNA_Output$RNA_log2FoldChange[k] <- round(TruncatedData$RNA_log2FoldChange[i], 3)
              RNA_Output$RNA_padj[k] <- round(TruncatedData$RNA_padj[i], 3)
              RNA_Output$RNA_pvalue[k] <- round(TruncatedData$RNA_pvalue[i],3)
          
            }}}
      }
      
    }
      included_patients <- unique(SubsetGenessplit$X1)
      
      
      ## Add the RPKM values for each gene fragment to the RNA_output
      fragment_RPKM_folder <- gsub('.{2}$', "", Truncated_DE_folder)
      fragment_RPKM_file <- list.files(fragment_RPKM_folder, full.names = F)[grep("RPKM", list.files(fragment_RPKM_folder, full.names = F))]
      
      fragment_RPKM <- read.delim(paste(fragment_RPKM_folder, fragment_RPKM_file, sep = ""), header = T, stringsAsFactors = F, check.names = F)
      # Some rownames contain erronous values, only select rows starting with "ENSG"
      fragment_RPKM <- fragment_RPKM[grep("ENSG", row.names(fragment_RPKM)),]
      # only select the truncated genes (wich have ensembl gene ids containing an underscore)
      fragment_RPKM <- fragment_RPKM[grep("_", row.names(fragment_RPKM)),]
      # 
      # if(grepl("MP", fragment_RPKM_file)){
      #   colnames(fragment_RPKM) <- gsub(".RNA.1", "", colnames(fragment_RPKM))
      #   colnames(fragment_RPKM) <- gsub("\\.", "-", colnames(fragment_RPKM))
      # }else if(grepl("TALK", fragment_RPKM_file)){
      #   # Remove the -RNA- suffices
      #   colnames(fragment_RPKM) <- gsub(".RNA.*", "", colnames(fragment_RPKM))
      #   
      #   # # The LCL cell lines have been sequenced in triplo, calculate the mean RPKM for each fragment 
      #   # fragment_RPKM <- t(fragment_RPKM)
      #   # fragment_RPKM <- aggregate(fragment_RPKM, by=list(rownames(fragment_RPKM)), mean)
      #   # rownames(fragment_RPKM) <- fragment_RPKM$Group.1
      #   # fragment_RPKM$Group.1 <- NULL 
      #   # fragment_RPKM <- data.frame(t(fragment_RPKM))
      #   # colnames(fragment_RPKM) <- gsub("X", "", colnames(fragment_RPKM))
      # }else{
      #   print("Unkown subset found")
      # }
      
      RNA_Output$ensembl_fragment <- ifelse(RNA_Output$SV_type == "Truncation",
                                            paste(RNA_Output$ensembl_gene_id, "_", RNA_Output$Fragment, sep = ""),
                                            RNA_Output$ensembl_gene_id)
      
      # Loop over each gene fragment
      for(i in 1:nrow(fragment_RPKM)){
        RPKMS_fragment <- fragment_RPKM[i,]
        
        # only add the RPKMs of the truncated gene fragment to the output if it is in RNA_output
        if(row.names(RPKMS_fragment) %in% RNA_Output$ensembl_fragment){
          data_to_replace <- RNA_Output[RNA_Output$ensembl_fragment == row.names(RPKMS_fragment),]
          #print(data_to_replace)
          # Some fragments occur in multiple patient, loop over each patient separately
          for(row in 1:nrow(data_to_replace)){
            # select the correct patient ID
            patient <- RNA_patients$RNA_files[RNA_patients$Patient_IDs == data_to_replace[row, "Patient"]]
            if(length(patient) > 0){
              rpkms_patient <- RPKMS_fragment[,c(grep(pattern = patient, x =  names(RPKMS_fragment)))]
              
              # Some samples have replicates (ie multiple columns in the RPKM file). Calculate the mean RPKM of the gene fragment for the samples
              if(!is.null(ncol(rpkms_patient))){
                RNA_Output$RNA_RPKM[RNA_Output$ensembl_fragment == row.names(RPKMS_fragment)][row] <- rowMeans(rpkms_patient)
              } else {
                RNA_Output$RNA_RPKM[RNA_Output$ensembl_fragment == row.names(RPKMS_fragment)][row] <- rpkms_patient
              }
              # Calculate the mean RPKM for all the other samples (= controls)
              RNA_Output$RNA_Controls_RPKM[RNA_Output$ensembl_fragment == row.names(RPKMS_fragment)][row] <- rowMeans(RPKMS_fragment[,-c(grep(pattern = patient, x =  names(RPKMS_fragment)))])
            }
            }

        }
      }
      
      # for(i in 1:nrow(RNA_Output)){
      #   patient <- RNA_patients$RNA_files[which(RNA_patients$Patient_IDs == RNA_Output$Patient[i])]
      #   #print(patient)
      #   if(length(patient) > 0){
      #     print(patient)
      #     if(patient %in% included_patients == TRUE ){
      #       ifelse(patient == "UTR22", patient <- "Pregno-child-blood", patient <- patient)
      #       gene <- RNA_Output$ensembl_fragment[i]
      #       col <- which(names(fragment_RPKM) == patient)
      #       row <- which(grepl(gene, row.names(fragment_RPKM)))
      #       rpkm <- fragment_RPKM[which(grepl(gene, row.names(fragment_RPKM))), which(names(fragment_RPKM) == patient)]
      #       cont_rpkm <- data.frame(rowMeans(fragment_RPKM[which(grepl(gene, row.names(fragment_RPKM))), which(names(fragment_RPKM) != patient)]))
      #       colnames(cont_rpkm) <- "RNA_Controls_RPKMs"
      #       ifelse(length(rpkm) == 0, next(), print(paste("add RPKM for:", gene, patient, sep = " ")))
      #       RNA_Output$RNA_RPKM[i] <- rpkm 
      #       RNA_Output$RNA_Controls_RPKM[i] <- cont_rpkm$RNA_Controls_RPKMs
      #     }else{
      #       next
      #     }}
      #   
  }
  return(RNA_Output)
}
