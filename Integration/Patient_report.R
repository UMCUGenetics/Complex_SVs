args <- commandArgs(trailingOnly=TRUE)
input_folder <- args[1]
ini_file <- args[2]
output_folder <- args[3]

#ini_file <- paste(input_folder, "/Scripts/Integration/Inputs/Genes_SVs.ini",sep = "")
# Read the ini file:
if(file.exists(ini_file)){
  ini <- read.delim(ini_file, header = T, stringsAsFactors = F, comment.char = "#")
  project_folder <- ini[which(ini$ID == "Project_folder"), "Value"]
  
} else {
  stop("Could not open ini file: ", ini_file)
}

library(reshape2)

source( paste(input_folder, "Scripts/Visualization/Breakpoint_browser.R", sep = ""))
source( paste(input_folder, "Scripts/Visualization/read_browser_data.R", sep = ""))
source( paste(input_folder, "Scripts/Visualization/Plot_ideograms_derivative_chromosomes.R", sep = ""))

#output_folder <- paste(input_folder, "SV_integrator/Results/MP/", sep = "")
Genes_SVs <- read.delim(paste(output_folder, "Genes_SVs.txt", sep = ""), stringsAsFactors = F)
standard.conf <- read.delim(paste(input_folder, "Visualization/conf/standard.conf", sep = ""), stringsAsFactors = F)

# Read the table containing the SVs:
SV_file <- ini[which(ini$ID == "SV_file"), "Value"]
if(file.exists(SV_file)){
  SVs <- read.delim(paste(project_folder,SV_file, sep = ""), stringsAsFactors = F)
} else {
  stop("Could not open SV_file: ", SV_file)
}

RNA_RPKM_files <- unlist(strsplit(ini[which(ini$ID == "RNA_RPKM_Files"), "Value"], split = ","))
RNA_DE_Folder <- ini[which(ini$ID == "RNA_DE_Folder"), "Value"]

if(!is.null(RNA_RPKM_files) & length(RNA_RPKM_files) > 0){
  RNA_counts <- list()
  for(RNA_counts_file in RNA_RPKM_files){
    
    print(paste("Reading: ", RNA_counts_file, sep= ""))
    RNA <- read.table(RNA_counts_file, check.names = F)
    
    RNA$ensembl_gene_ID <- row.names(RNA)
    RNA_counts[[RNA_counts_file]] <- RNA
  }
}

SV_heatmap <- function(Input_data,
                       Patient){
  
  Patient_data <- Input_data[which(Input_data$Patient == Patient),]
  
  # This number is used in the caption to show how many genes are not shown in the plot:
  Number_genes_patient <- length(unique(Patient_data$hgnc_symbol))
  
  # Filtering
  Patient_data <- Patient_data[which(Patient_data$phenoMatchScore > 0 |
                                       Patient_data$SV_type %in% c("Truncation", "Deletion", "Duplication") |
                                       !is.na(Patient_data$DDG2P) | Patient_data$pLI > 0.9 | Patient_data$RVIS < 0.1),]
  
  plot_data <- data.frame()
  
  # Facet 1 = gene info
  Patient_data$SV_type[Patient_data$SV_type == "Flanking_normal"] <- "Flanking_Normal"
  Patient_data$SV_type <- factor(Patient_data$SV_type)
  levels(Patient_data$SV_type) <- list( "DEL" = "Deletion", "DUP" = "Duplication", "INV" = "Inversion",
                                        "TRUNC" = "Truncation", "INTRA" = "Normal" , "INTRA" = "Translocation", "FLANK" = "Flanking_Inversion", "FLANK" = "Flanking_Normal",
                                        "FLANK" = "Flanking",
                                        "INS" = "Insertion")

  Patient_data$hgnc_symbol  <- ifelse(Patient_data$SV_type == "TRUNC", 
                                      paste(Patient_data$hgnc_symbol, Patient_data$Fragment, sep  = "_"), 
                                      Patient_data$hgnc_symbol )
  
  Patient_data$hgnc_symbol <- factor(Patient_data$hgnc_symbol,
                                     levels = Patient_data$hgnc_symbol[order(Patient_data$phenoMatchScore, decreasing = F, na.last = F)])
  
  ### Facet 1: gene info
  SV_Type <- data.frame(Gene = Patient_data$hgnc_symbol,
                        Variable = "SV Type",
                        Facet = "SV_info",
                        Value = ifelse(Patient_data$SV_type %in% c("TRUNC", "DEL", "DUP"), 50, -50),
                        Label = as.character(Patient_data$SV_type))
  plot_data <- rbind(plot_data, SV_Type)
  
  #plot_data$Gene <- factor(plot_data$Gene, levels = unique(as.vector(Patient_data$hgnc_symbol[order(Patient_data$phenoMatchScore, decreasing = F, na.last = F)])))
  
  Distance <- data.frame(Gene = Patient_data$hgnc_symbol,
                         Variable = "Distance\n(Mb)",
                         Facet = "SV_info",
                         Value = ifelse(Patient_data$Distance_to_BP < 1e6, (1e6 - Patient_data$Distance_to_BP) / 1e6 * 100 , 0),
                         Label = as.character(round(Patient_data$Distance_to_BP / 1e6, 2)))
  plot_data <- rbind(plot_data, Distance)

  
  ### Facet 2: Phenotype info
  
  Phenomatch_Score <- data.frame(Gene = Patient_data$hgnc_symbol,
                           Variable = "Phenomatch",
                           Facet = "Phenotype",
                           Value = ifelse(Patient_data$phenoMatchScore > 30, 100, Patient_data$phenoMatchScore/30*100),
                           Label = as.character(round(Patient_data$phenoMatchScore, 1)))
  plot_data <- rbind(plot_data, Phenomatch_Score)
  
  
  Phenomatch_Hits <- data.frame(Gene = Patient_data$hgnc_symbol,
                           Variable = "Phenohits",
                           Facet = "Phenotype",
                           Value = Patient_data$Phenomatches / Patient_data$HPO_Patient * 100,
                           Label = paste(Patient_data$Phenomatches, "/",  Patient_data$HPO_Patient, sep = ""))
  plot_data <- rbind(plot_data, Phenomatch_Hits)
  
  pLI <- data.frame(Gene = Patient_data$hgnc_symbol,
                                 Variable = "pLI",
                                 Facet = "Phenotype",
                                 Value = ifelse(Patient_data$pLI > 0.5, (Patient_data$pLI-0.5)*200 ,0),
                                 Label = as.character(round(Patient_data$pLI, 2)))
  plot_data <- rbind(plot_data, pLI)
  
  RVIS <- data.frame(Gene = Patient_data$hgnc_symbol,
                    Variable = "RVIS",
                    Facet = "Phenotype",
                    Value = ifelse(Patient_data$RVIS < 50, 100 - Patient_data$RVIS*2,0),
                    Label = as.character(round(Patient_data$RVIS, 1)))
  plot_data <- rbind(plot_data, RVIS)
  
  

  
  Patient_data$DDG2P <- factor(Patient_data$DDG2P)
  levels(Patient_data$DDG2P) <- list("Yes" = "confirmed", "Prob" = "probable", "Pos"= "possible", "No" = NA)
  DDG2P <- data.frame(Gene = Patient_data$hgnc_symbol,
                     Variable = "DDG2P",
                     Facet = "Phenotype",
                     Value = ifelse(Patient_data$DDG2P %in% c("Yes", "Prob"), 75, 0),
                     Label = as.character(Patient_data$DDG2P))
  plot_data <- rbind(plot_data, DDG2P)
  
  HI <- data.frame(Gene = Patient_data$hgnc_symbol,
                   Variable = "HI",
                   Facet = "Phenotype",
                   Value = ifelse(Patient_data$HI < 50, 100 - Patient_data$HI*2, 0),
                   Label = as.character(round(Patient_data$HI, 1)))
  plot_data <- rbind(plot_data, HI)
  
  inheritance <- data.frame(Gene = Patient_data$hgnc_symbol,
                   Variable = "Inheritance",
                   Facet = "Phenotype",
                   Value = ifelse(Patient_data$Inheritance %in% c("AD", "AD+AR"), 75, 0),
                   Label = as.character(Patient_data$Inheritance))
  plot_data <- rbind(plot_data, inheritance)

  ### Facet 3: RNA+ATAC
  if(!is.null(Patient_data$RNA_RPKM)){
    RNA <- data.frame(Gene = Patient_data$hgnc_symbol,
                      Variable = "RNA",
                      Facet = "Expression",
                      Value = ifelse(Patient_data$RNA_RPKM > 0.5 | Patient_data$RNA_Controls_RPKM > 0.5, 
                                     Patient_data$RNA_log2FoldChange / 2 * -100,
                                     NA),
                      Label = ifelse(Patient_data$RNA_RPKM > 0.5 | Patient_data$RNA_Controls_RPKM > 0.5,
                                     ifelse(Patient_data$RNA_pvalue <0.05, paste(as.character(round(Patient_data$RNA_log2FoldChange, 1)), "*",sep= ""),
                                            as.character(round(Patient_data$RNA_log2FoldChange, 1))),
                                     "NE"))
    RNA$Value <- ifelse(RNA$Value > 100, 100, RNA$Value)
    RNA$Value <- ifelse(RNA$Value < -100, -100, RNA$Value)
    plot_data <- rbind(plot_data, RNA)
  }

  
  ### Facet 4: 3D genome
  TADs <- data.frame(Gene = Patient_data$hgnc_symbol,
                     Variable = "TADs",
                     Facet = "Disrupted interactions",
                     Value = Patient_data$Average_TAD_Loss,
                     Label = as.character(Patient_data$Cell_Types_TAD_Disrupted))
  
  if(nrow(TADs) > 0){
    plot_data <- rbind(plot_data, TADs)
  }
  
  V4C <- data.frame(Gene = Patient_data$hgnc_symbol,
                     Variable = "Virtual 4C",
                     Facet = "Disrupted interactions",
                     Value = Patient_data$V4C_Disruption,
                     Label = as.character(round(Patient_data$V4C_Disruption, 0)))
  
  if(nrow(V4C) > 0){
    plot_data <- rbind(plot_data, V4C)
  }
  
  PCHiC <- data.frame(Gene = Patient_data$hgnc_symbol,
                      Variable = "PCHiC",
                      Facet = "Disrupted interactions",
                      Value = Patient_data$PCHiC_Disruptions / Patient_data$PCHiC_Interactions * 100,
                      Label = paste(Patient_data$PCHiC_Disruptions, "/", Patient_data$PCHiC_Interactions,sep = ""))
  if(nrow(PCHiC) > 0){
    plot_data <- rbind(plot_data, PCHiC)
  }
  
  DHS <- data.frame(Gene = Patient_data$hgnc_symbol,
                      Variable = "DHS",
                      Facet = "Disrupted interactions",
                      Value = Patient_data$Disrupted_DHS / Patient_data$Distal_DHS_connections * 100,
                      Label = paste(Patient_data$Disrupted_DHS, "/", Patient_data$Distal_DHS_connections,sep = ""))
  if(nrow(DHS) > 0){
    plot_data <- rbind(plot_data, DHS)
  }
  
  Enhancers <- data.frame(Gene = Patient_data$hgnc_symbol,
                    Variable = "Enhancers",
                    Facet = "Disrupted interactions",
                    Value = Patient_data$Enhancer_loss,
                    Label = round(Patient_data$Enhancer_loss, 1), stringsAsFactors = F)
  if(nrow(Enhancers) > 0){
    plot_data$Label <- as.character(plot_data$Label)
    plot_data <- rbind(plot_data, Enhancers)
  }
  
  Final_score <- data.frame(Gene = Patient_data$hgnc_symbol,
                          Variable = "Score",
                          Facet = "Disrupted interactions",
                          Value = ifelse(Patient_data$Score < 100, Patient_data$Score, 100),
                          Label = round(Patient_data$Score, 1), stringsAsFactors = F)
  if(nrow(Final_score) > 0){
      plot_data <- rbind(plot_data, Final_score)
  }
  
  

  # facet 5: classication
  if(is.null(Patient_data$Conclusion) == FALSE){
    Patho <- data.frame(Gene = Patient_data$hgnc_symbol,
                        Variable = "Driver",
                        Facet = "Class",
                        Value = ifelse(Patient_data$Conclusion == "Strong_driver", 50, 0),
                        Label = gsub(pattern = "_driver", replacement = "", x = Patient_data$Conclusion))
    Patho$Value <- ifelse(Patho$Label == "Weak", 25,  Patho$Value)
    #levels(Patho$Label) <- list("LP" = "Likely_pathogenic", "LB" = "Likely_benign", "LBA" = "Likely_benign_affected", "UP" = "Unlikely_pathogenic", "PP" = "Possibly_pathogenic")
    plot_data <- rbind(plot_data, Patho)
  }
  
  if(is.null(Patient_data$manual_classification) == FALSE){
    Patho_manual <- data.frame(Gene = Patient_data$hgnc_symbol,
                        Variable = "Classification",
                        Facet = "Class",
                        Value = ifelse(Patient_data$manual_classification == "Likely_pathogenic", 50, 0),
                        Label = Patient_data$manual_classification)
    Patho_manual$Value <- ifelse(Patho_manual$Label == "Possibly_pathogenic", 25,  Patho_manual$Value)
    levels(Patho_manual$Label) <- list("LP" = "Likely_pathogenic", "LB" = "Likely_benign", "LBA" = "Likely_benign_affected", "UP" = "Unlikely_pathogenic", "PP" = "Possibly_pathogenic")
    plot_data <- rbind(plot_data, Patho_manual)
  }
  
  plot_data$plot <- 1
  
  
  plot_data2 <- plot_data
  
  if(length(unique(plot_data$Gene)) > 60){
    
    plot_data$plot[which(plot_data$Gene %in% levels(plot_data$Gene)[1:(floor(length(levels(plot_data$Gene))/2)-1)])] <- 2
    
    #plot_data1 <- plot_data[which(plot_data$Gene %in% levels(plot_data$Gene)[floor(length(levels(plot_data$Gene))/2):(length(levels(plot_data$Gene)))]),]
    #plot_data2 <- plot_data[which(plot_data$Gene %in% levels(plot_data$Gene)[1:(floor(length(levels(plot_data$Gene))/2)-1)]),]
    
  }
  
  caption <- paste("Showing data for ", length(unique(Patient_data$hgnc_symbol)), "/", Number_genes_patient, " genes +/- 2 Mb from SVs", sep = "")
  
  bottommargin <- ifelse(length(unique(plot_data$Gene)) < 50, 23 - length(unique(plot_data$Gene)) * 0.5, 0.5)
  
  for(plot_nr in unique(plot_data$plot)){
    print(plot_nr)
    plot <- ggplot(plot_data[which(plot_data$plot == plot_nr),], aes(x = Variable, y = Gene, fill = Value)) + 
      geom_tile(col = "black") + 
      geom_text(aes(label = Label), size = 3) + 
      facet_grid(~ Facet, scales = "free", space = "free") + 
      scale_fill_gradient2(low = "darkblue", mid = "white", high = rgb(178,34,34, maxColorValue = 255), na.value = "gray")+
      scale_x_discrete(position = "top") +
      theme_minimal()+
      theme(axis.ticks.y=element_blank(), 
            axis.title.x=element_blank(),
            axis.text.x = element_text(angle = 60, hjust = 0, size = 12), plot.title = element_text(hjust = 0.5),
            strip.placement = "outside", strip.background = element_rect(fill = "lightgray", colour = "white"), strip.text = element_text(face = "bold", size = 10), 
            plot.margin = unit(c(0.5,0.5,bottommargin,0.5), "cm")) + 
      labs(caption = caption) + ggtitle(label = ifelse(plot_nr != 2, Patient, paste(Patient, " (2) ", sep = ""))) 
    
    print(plot)
    
  }

}



plot_phenomatch <- function(Patient,
                            phenomatch_folder =  paste(input_folder, "Results/MP/Phenomatch/", sep = ""),
                            Genes_SVs = Gene_data,
                            SVs = SVs){
  
  phenomatch_file_patient <- paste(phenomatch_folder, "Phenomatch_", Patient, ".txt", sep = "")
  
  phenomatch_overview <- read.delim(phenomatch_file_patient, check.names = F, stringsAsFactors = F)
  
  patient_data <- Genes_SVs[which(Genes_SVs$Patient == Patient),]
  patient_data$gene_symbol <- patient_data$hgnc_symbol
  
  # Only select one row per gene (removes one of the truncated gene fragments)
  patient_data <- patient_data[!duplicated(patient_data$gene_symbol),]
  
  phenomatch_overview <- merge(phenomatch_overview, patient_data[,c("hgnc_symbol", "SV_type")], by = "hgnc_symbol", all.y = F)
  
  phenomatch_overview$SV_type <- factor(phenomatch_overview$SV_type)
  levels(phenomatch_overview$SV_type) <- list( "DEL" = "Deletion", "DUP" = "Duplication", "INV" = "Inversion",
                                               "TRUNC" = "Truncation", "INTRA" = "Normal" , "FLANK" = "Flanking_Inversion", "FLANK" = "Flanking_Normal",
                                               "INS" = "Insertion")
  
  # Phenomatch automatically looks at genes +/- 3Mb from the SVs and the output therefore contains genes we're not looking at.
  # Select the genes +/- 2Mb from the SVs
  
  phenomatch_overview_2Mb <- phenomatch_overview[which(phenomatch_overview$hgnc_symbol %in% as.vector(patient_data$hgnc_symbol)),]
  phenomatch_overview_2Mb$hgnc_symbol <- factor(phenomatch_overview_2Mb$hgnc_symbol, 
                                                levels = phenomatch_overview_2Mb$hgnc_symbol[order(phenomatch_overview_2Mb$phenoMatchScore, decreasing = F)])
  
  phenomatch_overview_melted <- melt(phenomatch_overview_2Mb, id.vars = "hgnc_symbol")
  #phenomatch_overview_melted$value <- as.numeric(phenomatch_overview_melted$value)
  phenomatch_overview_melted$Score <- phenomatch_overview_melted$value
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "SV_type"] <- 0
  phenomatch_overview_melted$Score[phenomatch_overview_melted$value %in% c("DEL", "DUP", "TRUNC")] <- 15
  
  phenomatch_overview_melted$Score <- as.numeric(phenomatch_overview_melted$Score)
  phenomatch_overview_melted$Score <- ifelse(phenomatch_overview_melted$Score > 30, 30, phenomatch_overview_melted$Score)
  
  phenomatch_overview_melted$grid <- ifelse(phenomatch_overview_melted$variable == "SV_type", "SV", "Phenomatch")
  
  phenomatch_overview_melted$grid <- factor(phenomatch_overview_melted$grid, levels = c("SV", "Phenomatch"))
  
  phenomatch_overview_melted$value[phenomatch_overview_melted$variable != "SV_type"] <- round(as.numeric(phenomatch_overview_melted$value[phenomatch_overview_melted$variable != "SV_type"]), 1)
  
  names(phenomatch_overview_melted) <- c("Gene", "Phenotype", "phenoMatchScore", "Score", "grid")
  
  output_height <- nrow(phenomatch_overview) / 3 + 2
  
  #print(output_height)
  output_width <- ncol(phenomatch_overview) / 2 + 2
  #print(output_width)
  
  caption <- paste("Showing data for ", length(unique(phenomatch_overview_melted$Gene)), "/", length(unique(patient_data$hgnc_symbol)), " genes +/- 2 Mb from SVs", sep = "")
  
  bottommargin <-  ifelse(length(unique(phenomatch_overview_melted$Gene)) < 50, 20 - length(unique(phenomatch_overview_melted$Gene)) * 0.5, 0.5)
  
  data_table <- ggplot(phenomatch_overview_melted, aes(x = Phenotype, y = Gene, fill = Score, label = phenoMatchScore)) + 
    geom_tile(lwd = 0.4, colour  = "black") +
    geom_text(size = 3, data = phenomatch_overview_melted, aes(label = phenoMatchScore)) + 
    scale_fill_gradient(low ="white", high = "red", limits=c(0,30), na.value = "gray") +
    scale_x_discrete(position = "top") + 
    theme_minimal(base_size = 10) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 0), axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5),
          strip.placement = "outside", strip.background = element_rect(fill = "lightgray", colour = "white"), strip.text = element_text(face = "bold", size = 10),
          plot.margin = unit(c(0.5,0.5,bottommargin,0.5), "cm")) + 
    ggtitle(label = Patient) +
    facet_grid(~grid, scales = "free", space = "free") +
    labs(caption = caption) 
  
  print(data_table)

}




plot_RNA_boxplots <- function(patient,
                              RNA_data = RNA_counts,
                              Genes_SVs,
                              max_number_genes = 20){
  
  
  Genes_SVs_Patient <- Genes_SVs[which(Genes_SVs$Patient == patient),]
  
  Genes_SVs_Patient <- Genes_SVs_Patient[which(Genes_SVs_Patient$RNA_RPKM > 0.5 | Genes_SVs_Patient$RNA_Controls_RPKM > 0.5),]
  Genes_SVs_Patient <- Genes_SVs_Patient[which(Genes_SVs_Patient$phenoMatchScore > 0 | Genes_SVs_Patient$RNA_pvalue < 0.05),]
  
  
  if(nrow(Genes_SVs_Patient) > max_number_genes){
    Genes_SVs_Patient <- Genes_SVs_Patient[1:max_number_genes,]
  }
  
  for(counts in names(RNA_data)){

    patients_counts <- gsub("-RNA-[[:digit:]]$", "", names(RNA_data[[counts]]))
    patients_counts <- gsub("MP[0-9][0-9]-", "MP", patients_counts)
    
    patients_counts <- ifelse(grepl(patients_counts, pattern = "Pregno-child") == TRUE, "UTR22", patients_counts)

    
    if(patient %in% patients_counts){
      
      RNA_counts_patient <- RNA_data[[counts]]
    }
  }
  
  Genes_to_plot <- Genes_SVs_Patient[,c("ensembl_gene_id",  "hgnc_symbol", "RNA_pvalue", "SV_type")]
  Genes_to_plot$Label <- ifelse(Genes_to_plot$RNA_pvalue < 0.05, paste(Genes_to_plot$hgnc_symbol, "*", sep = ""), Genes_to_plot$hgnc_symbol)
  Genes_to_plot$Label <- ifelse(Genes_to_plot$SV_type == "Deletion", paste(Genes_to_plot$Label, " (d)", sep = ""), Genes_to_plot$Label)
  Genes_to_plot$Label <- ifelse(Genes_to_plot$SV_type == "Duplication", paste(Genes_to_plot$Label, " (D)", sep = ""), Genes_to_plot$Label)
  
  names(Genes_to_plot)[names(Genes_to_plot) == "ensembl_gene_id"] <- "ensembl_gene_ID"
  Genes_to_plot <- merge(Genes_to_plot, RNA_counts_patient, by= "ensembl_gene_ID")
  Genes_to_plot <- Genes_to_plot[,-which(names(Genes_to_plot) %in% c("ensembl_gene_ID", "SV_Type","hgnc_symbol", "RNA_pvalue"))]
  
  
  plot_data <- melt(Genes_to_plot)
  plot_data$variable <- gsub("-RNA-[[:digit:]]$", "", plot_data$variable)
  plot_data$variable <- gsub("MP[0-9][0-9]-", "MP", plot_data$variable)
  plot_data$variable <- ifelse(grepl(plot_data$variable, pattern = "Pregno-child") == TRUE, "UTR22", plot_data$variable)
  
  plot_data$type <- ifelse(grepl(pattern = patient, plot_data$variable), patient, "Control")
  plot_data$type <- factor(plot_data$type, levels = c("Control", patient))
  
  number_rows <- ceiling(length(unique(plot_data$Label)) / 4)
  
  bottommargin <-  ifelse(number_rows < 5, (5-number_rows)*5, 0.5)
  
  if(nrow(plot_data) > 0){
    plot <- ggplot(plot_data, aes(x = "Blood", y = value)) +
      stat_boxplot(geom = "errorbar", width = 0.09, lwd = 0.65, col = "gray40") +
      geom_boxplot(col = "gray50", fill = "gray90", outlier.shape = NA, lwd = 0.6) + 
      geom_jitter(aes(x = "Blood", y= value, col = type, size = type), width= 0.25, alpha = 0.9) + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5),
            plot.margin = unit(c(0.5,0.5,bottommargin,0.5), "cm"),
            panel.grid.major.x = element_blank()) +
      facet_wrap(~ Label, scales = "free", ncol = 4) +
      scale_color_manual(values = c("dodgerblue","darkred"))  + 
      expand_limits(y = 0) +
      scale_size_discrete(range=c(1,3)) +
      labs(x = "", y = "RNA (RPKM)") + ggtitle(paste("RNA expression ", patient, sep = ""))
    
    print(plot)
  }
  
}


plot_truncated_gene <- function(ensembl_gene_id, SVs, patient, exon_height = 0.6, dis_between_transc = 0.2){
  ensembl_hg19 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  
  refseqs_gene <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id", "ensembl_transcript_id","refseq_mrna"),
                        filters = c("ensembl_gene_id"), values = ensembl_gene_id, mart = ensembl_hg19)
  
  refseqs_gene <- refseqs_gene[which(refseqs_gene$refseq_mrna != ""),]
  
  exons_gene <- getBM(attributes = c( "ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand"),
                      filters = c("ensembl_transcript_id"), values = refseqs_gene$ensembl_transcript_id, mart = ensembl_hg19)
  
  exons_gene <- merge(exons_gene, refseqs_gene[,c("ensembl_transcript_id", "refseq_mrna", "hgnc_symbol")])
  print(head(exons_gene))
  xlim_plot <- c(min(exons_gene$exon_chrom_start)-0.5e5, max(exons_gene$exon_chrom_end)+0.5e5)
  
  height_plot <- length(unique(exons_gene$refseq_mrna)) * (exon_height+dis_between_transc) + 1
  
  ylim <- c(0,height_plot)
  
  par(mar=c(1, 6, 1, 1))
  #par(mar=c(1, 5, 2, 1))
  plot.new()
  
  plot.window(xlim=xlim_plot, ylim=ylim)
  
  y0 = 0.5
  
  for(transcript in unique(exons_gene$refseq_mrna)){
    print(transcript)
    
    exons_transcripts <- exons_gene[which(exons_gene$refseq_mrna == transcript),]
    print(head(exons_transcripts))
    
    if(max(exons_transcripts$exon_chrom_start) - min(exons_transcripts$exon_chrom_start) > 20000){
      arrows(x0 = seq(min(exons_transcripts$exon_chrom_start)+10000, max(exons_transcripts$exon_chrom_start) - 10000, by = 10000),
             x1 = seq(min(exons_transcripts$exon_chrom_start)+10000, max(exons_transcripts$exon_chrom_start) - 10000, by = 10000)+1000,
             y0 = y0+exon_height/2, code = ifelse(unique(exons_transcripts$strand == 1), 2, 1),
             col = "darkgray", length = 0.08, angle = 60)
    }
    
    arrows(x0 = min(exons_transcripts$exon_chrom_start), x1 = max(exons_transcripts$exon_chrom_start), y0 = y0+exon_height/2, code = 0, lwd = 1)
    
    rect(xleft = exons_transcripts$exon_chrom_start, xright = exons_transcripts$exon_chrom_end, ybottom = y0, ytop = y0+exon_height, col = "darkgray", lwd  = 1)
    
    mtext(text = transcript, side = 2, at = y0+exon_height/2, las = 1, line = -2, cex = 0.8)
    y0  <- y0 + exon_height + dis_between_transc
    
  }
  
  SVs_patient <- SVs[which(SVs$Patient == patient & SVs$chr == unique(exons_gene$chromosome_name)),]
  arrows(x0 = SVs_patient$end, y0= 0, y1 = 10, col = rgb(128, 0,0,150, maxColorValue = 255), lty = 2, code = 0, lwd = 1.5)
  print(y0)
  text(x = SVs_patient$end - 0.2e5, y = y0+exon_height/2, labels = SVs_patient$Fragment, col = "red")
  text(x = SVs_patient$start + 0.2e5, y = y0+exon_height/2, labels = SVs_patient$Fragment, col = "red")
  
  
  axis(1, line = -1)
  title(paste(unique(exons_gene$hgnc_symbol), " (", patient, ")", sep = ""))
}



dir.create(output_folder, showWarnings = FALSE)

for(patient in unique(Genes_SVs$Patient)){
  if(patient != "17078"){
    print(paste("# Generating a report for patient: ", patient, sep = ""))
    
    pdf(file = paste(output_folder, "Report_", patient, "-2.pdf", sep = ""), width = 8, height =  11, pointsize = 10)
    
    SVs_patient <- SVs[which(SVs$Patient == patient),]
    
    if(nrow(SVs_patient) > 25){
      par(mfrow = c(1,1))
    } else {
      par(mfrow = c(2,1))
    }
    plot_der_karyotype(Patient = patient, breakpoints = SVs)
    box("figure")  
    plot_der_chroms(Patient = patient, breakpoints=SVs)
    box("figure")  
    par(mfrow = c(1,1))
    
    SV_heatmap(Input_data = Genes_SVs, Patient = patient)
    
    plot_phenomatch(Genes_SVs = Genes_SVs, Patient = patient, phenomatch_folder =  paste(input_folder, "SV_integrator/Results/MP/Phenomatch/", sep = ""))
    #par(mfrow = c(1,1), mar = c(10, 1,1,1))
    
    plot_regions <- define_plot_region(patient = patient, derivative = FALSE, SVs = SVs)
    panels <- length(plot_regions)/3
    par(mfrow = c(2,1))
    for(panel in 0:(panels-1)){
      print(panel)
      xlim <- c(as.numeric(plot_regions[2 + panel*3]),as.numeric(plot_regions[3 + panel*3]))
      print(xlim)
      
      breakpoint_plot(conf = standard.conf, patient = patient, chromosome = plot_regions[1 + panel*3],
                    xlim_plot = c(as.numeric(plot_regions[2 + panel*3]),as.numeric(plot_regions[3 + panel*3])), Genes_SVs = Genes_SVs)
      box("figure") 
    }
    
    # plot_der_regions <- define_plot_region(patient = patient, derivative = TRUE, SVs = SVs)
    # panels_der <- length(plot_der_regions)/3
    # for(panel in 0:(panels_der-1)){
    #   print(panel)
    #   xlim <- c(as.numeric(plot_der_regions[2 + panel*3]),as.numeric(plot_der_regions[3 + panel*3]))
    #   print(xlim)
    #   
    #   breakpoint_plot(conf = standard.conf, patient = patient, Genes_SVs = Genes_SVs, chromosome = plot_der_regions[1 + panel*3],
    #                   xlim_plot = c(as.numeric(plot_der_regions[2 + panel*3]),as.numeric(plot_der_regions[3 + panel*3])))
    #   box("figure") 
    # }
    par(mfrow = c(1,1))
    
    if(!is.na(mean(Genes_SVs[Genes_SVs$Patient == patient, "RNA_RPKM"]))){
      print("Plotting RNA-seq data")
      plot_RNA_boxplots(patient = patient, Genes_SVs = Genes_SVs)
    } else {
      print(paste("No RNA data for patient ", patient, sep =""))
    }
  
    
    par(mfrow = c(4,2))
    
    for(truncated_gene in unique(Genes_SVs$ensembl_gene_id[Genes_SVs$SV_Type == "Truncation" & Genes_SVs$Patient == patient])){
      plot_truncated_gene(ensembl_gene_id = truncated_gene, patient = patient, SVs = SVs)
      box(which = "figure")
    }

    dev.off()
  }

}

#Genelist <- read.delim(file = ini[which(ini$ID == "Genelist"), "Value"])
#Exons <- read.delim(file = ini[which(ini$ID == "Exons"), "Value"])


PCHiC_Disruption_plot <- function(PCHiC_overview = test){
  PCHiC_overview$Patient <- sub("_.*", "", PCHiC_overview$Gene_ID)
  
  for(patient in unique(PCHiC_overview$Patient)){
    print(patient)
    plot_data <- data.frame()
    
    patient_data <- PCHiC_overview[which(PCHiC_overview$Patient == patient),]
    
    
    for(cell_type in PCHiC_Celltypes){
      print(cell_type)
      Patient_cell_data <- data.frame(Gene = gsub(".*_", "", patient_data$Gene_ID))
      Patient_cell_data$Cell_type <- cell_type
      Patient_cell_data$Tissue <- Cells_metadata[which(Cells_metadata$Cell_type == cell_type), "Tissue"]
      
      losses <- patient_data[,grep(paste("\\b", cell_type, "_Disrupted", sep = ""), x = names(patient_data))]
      interactions <- patient_data[,grep(paste("\\b", cell_type, "_Interactions", sep = ""), x = names(patient_data))]
      Patient_cell_data$Interactions <- interactions
      Patient_cell_data$PCHiC_Disruptions <- (losses / interactions) * 100
      plot_data <- rbind(plot_data, Patient_cell_data)
      
    }
    
    ggplot(plot_data, aes(x = Cell_type, y = Gene, fill = PCHiC_Disruptions)) + 
      geom_tile() + geom_text(aes(label = Interactions)) + 
      facet_grid(~Tissue, scales = "free", space = "free") + 
      scale_fill_gradient(low = "white", high = "darkred", na.value = "lightgray") + ggtitle(paste("PCHiC Disruptions ", patient, sep = ""))+
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  }
}
